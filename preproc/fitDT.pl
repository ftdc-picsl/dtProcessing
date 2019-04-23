#!/usr/bin/perl -w
#
# Fit DTI to corrected DWI data
#


use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;
use File::Basename;
use Getopt::Long;


my $usage = qq{

  $0
     --dwi 
     --bvecs
     --bvals
     --mask 
     --output-root 
     [options]

  
  Fit the diffusion tensor to some data.


  Required args:

   --dwi
     4D DWI data, from which the b=0 data is extracted.
  
   --bvecs
     bvecs for the data.

   --bvals
     bvals for the data,

   --mask
     Brain mask, tensors fitted only within the mask.

   --output-root
     Output root, prepended to output file names.
  
  Options:

   --restore-sigma
     Outlier rejection sigma, only used for RESTORE. Larger values mean fewer measurements rejected as outliers. The
     noise standard deviation, estimated from the b=0 data, is a reasonable lower bound for sigma.

   --algorithm 
     Fitting algorithm for the tensor, choice of "dt" (linear least squares), "wdt" (weighted linear least squares), or "RESTORE" (default = wdt)


  Requires ANTs and Camino

};


if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

# Get the directories containing programs we need
my ($antsPath, $sysTmpDir) = @ENV{'ANTSPATH', 'TMPDIR'};


if (!$antsPath || ! -f "${antsPath}antsRegistration") {
    die("Script requires ANTSPATH to be defined");
}

my $dtfit = `which dtfit`;

chomp $dtfit;

if (! -f $dtfit) {
    die("Script requires Camino");
}


my ($dwiData, $bvecs, $bvals, $mask, $outputRoot);

my $restoreSigma = 0;
my $algorithm = "wdt";


GetOptions ("dwi=s" => \$dwiData,  
            "bvals=s" => \$bvals,
            "bvecs=s" => \$bvecs,
	    "mask=s" => \$mask,
            "output-root=s" => \$outputRoot,
	    "algorithm=s" => \$algorithm,
	    "restore-sigma" => \$restoreSigma
    )
    or die("Error in command line arguments\n");


$algorithm = lc($algorithm);

# Treat b-values less than this as zero
my $effectiveB0 = 10;

my ($outputFileRoot,$outputDir) = fileparse($outputRoot);

if (! -d $outputDir ) { 
  mkpath($outputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputDir\n\t";
}

# Set to 1 to delete intermediate files after we're done
# Has no effect if using qsub since files get cleaned up anyhow
my $cleanup=1;

# Directory for temporary files that is deleted later if $cleanup
my $tmpDir = "";

my $tmpDirBaseName = "${outputFileRoot}DTIfit";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = $outputDir . "/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = $sysTmpDir . "/${tmpDirBaseName}";
}

# Gets removed later, so check we can create this and if not, exit immediately
mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir (maybe it exists from a previous failed run)\n\t";

my $scheme = "${outputRoot}MergedSchemeCorrected.scheme";

# Make scheme file
system("fsl2scheme -bvals $bvals -bvecs $bvecs -bscale 1 -outputfile $scheme");

# Get number of measurements from bvals, used to reshape output
my $numMeas=`wc -w $bvals | cut -d ' ' -f 1`;

chomp $numMeas;

# output some diagnostic stuff

# Average B0 - used to compute SNR image if we have sigma from weighted fit
system("averagedwi -inputfile $dwiData -outputfile ${outputRoot}AverageB0.nii.gz -minbval 0 -maxbval $effectiveB0 -schemefile $scheme");

# Average DWI is mostly for visual evaluation - output average of the shell in the vicinity of 500-1600
system("averagedwi -inputfile $dwiData -outputfile ${outputRoot}AverageDWI.nii.gz -minbval 500 -maxbval 1600 -schemefile $scheme");


if ($algorithm eq "dt") {
    fitDT($dwiData, $scheme, $mask, $outputRoot);
}
elsif ($algorithm eq "wdt") {
    fitWeightedDT($dwiData, $scheme, $mask, $outputRoot);
}
elsif ($algorithm eq "restore") {
    fitRestoreDT($dwiData, $scheme, $mask, $restoreSigma, $outputRoot);
}
else {
    die("Unrecognized tensor fitting algorithm $algorithm");
}


# Now compute FA, AD, RD, MD etc with ANTs
system("${antsPath}ImageMath 3 ${outputRoot}FA.nii.gz TensorFA ${outputRoot}DT.nii.gz");
system("${antsPath}ImageMath 3 ${outputRoot}MD.nii.gz TensorMeanDiffusion ${outputRoot}DT.nii.gz");
system("${antsPath}ImageMath 3 ${outputRoot}RD.nii.gz TensorRadialDiffusion ${outputRoot}DT.nii.gz");
system("${antsPath}ImageMath 3 ${outputRoot}AD.nii.gz TensorAxialDiffusion ${outputRoot}DT.nii.gz");
system("${antsPath}ImageMath 3 ${outputRoot}RGB.nii.gz TensorColor ${outputRoot}DT.nii.gz");


# cleanup

if ($cleanup) {
    system("rm ${tmpDir}/*");
    system("rmdir ${tmpDir}");
} 


# Output in a consistent style, not Camino underscore convention. Eg, ${outputRoot}DT.nii.gz

sub fitDT {

    my ($data, $scheme, $mask, $outputRoot) = @_; 

    system("dtfit $data $scheme -brainmask $mask > ${tmpDir}/dt.Bdouble");

    system("dtspd -inputfile ${tmpDir}/dt.Bdouble -dwifile $data -schemefile $scheme -brainmask $mask -unweightedb $effectiveB0 -editeigenvalues -outputroot ${tmpDir}/spd_");

    system("cp ${tmpDir}/spd_brainMask.nii.gz ${outputRoot}BrainMask.nii.gz");
    system("cp ${tmpDir}/spd_numNegDiffCoeff.nii.gz ${outputRoot}NumNegDiffCoeff.nii.gz");
    system("cp ${tmpDir}/spd_numNegEV.nii.gz ${outputRoot}NumNegEV.nii.gz");

    system("dt2nii -inputfile ${tmpDir}/spd_dtSPD.Bdouble -outputdatatype float -gzip -header $data -outputroot ${tmpDir}/spdNii_");

    system("cp ${tmpDir}/spdNii_exitcode.nii.gz ${outputRoot}ExitCode.nii.gz");
    system("cp ${tmpDir}/spdNii_dt.nii.gz ${outputRoot}DT.nii.gz");

}


sub fitWeightedDT {

    my ($data, $scheme, $mask, $outputRoot) = @_; 

    # sigmaSq is always written as double, the tensor type is controlled by -outputdatatype
    # Add -residualmap ${tmpDir}/residuals.Bdouble if you want residuals, but this is a lot of disk space
    #
    system("wdtfit $data $scheme ${tmpDir}/sigmaSq.Bdouble -brainmask $mask > ${tmpDir}/dt.Bdouble");
    
    system("cat ${tmpDir}/sigmaSq.Bdouble | voxel2image -inputdatatype double -outputdatatype float -header $data -outputvector -outputroot ${tmpDir}/SigmaSq");

    system("${antsPath}ImageMath 3 ${tmpDir}/Sigma.nii.gz ^ ${tmpDir}/SigmaSq.nii.gz 0.5");

    # Sigma probably more useful for QC purposes
    system("cp ${tmpDir}/Sigma.nii.gz ${outputRoot}Sigma.nii.gz");

    # If we have residuals, convert to nii
    # system("cat ${tmpDir}/residuals.Bdouble | voxel2image -inputdatatype double -outputdatatype float -header $data -outputvector -components $numMeas -outputroot ${outputRoot}StandardizedResiduals");
    system("dtspd -inputfile ${tmpDir}/dt.Bdouble -dwifile $data -schemefile $scheme -brainmask $mask -editeigenvalues -outputroot ${tmpDir}/spd_");

    system("cp ${tmpDir}/spd_brainMask.nii.gz ${outputRoot}BrainMask.nii.gz");
    system("cp ${tmpDir}/spd_numNegDiffCoeff.nii.gz ${outputRoot}NumNegDiffCoeff.nii.gz");
    system("cp ${tmpDir}/spd_numNegEV.nii.gz ${outputRoot}NumNegEV.nii.gz");

    system("dt2nii -inputfile ${tmpDir}/spd_dtSPD.Bdouble -outputdatatype float -gzip -header $data -outputroot ${tmpDir}/spdNii_");

    system("cp ${tmpDir}/spdNii_exitcode.nii.gz ${outputRoot}ExitCode.nii.gz");
    system("cp ${tmpDir}/spdNii_dt.nii.gz ${outputRoot}DT.nii.gz");

}


sub fitRestoreDT {

    my ($data, $scheme, $mask, $restoreSigma, $outputRoot) = @_; 

    system("restore $data $scheme $restoreSigma ${tmpDir}/outliers.Bbyte -outputdatatype float -brainmask $mask > ${tmpDir}/dt.Bfloat");
    
    system("cat ${tmpDir}/outliers.Bbyte | voxel2image -inputdatatype byte -header $data -outputvector -components $numMeas -outputroot ${outputRoot}RestoreOutliers");

    # Keep original exit code as it contains outlier information
    system("cp ${tmpDir}/dt_exitcode.nii.gz ${outputRoot}RestoreExitCode.nii.gz");    

    system("dtspd -inputfile ${tmpDir}/dt.Bdouble -dwifile $data -schemefile $scheme -brainmask $mask -editeigenvalues -outputroot ${tmpDir}/spd_");

    system("cp ${tmpDir}/spd_brainMask.nii.gz ${outputRoot}BrainMask.nii.gz");
    system("cp ${tmpDir}/spd_numNegDiffCoeff.nii.gz ${outputRoot}NumNegDiffCoeff.nii.gz");
    system("cp ${tmpDir}/spd_numNegEV.nii.gz ${outputRoot}NumNegEV.nii.gz");

    system("dt2nii -inputfile ${tmpDir}/spd_dtSPD.Bdouble -outputdatatype float -gzip -header $data -outputroot ${tmpDir}/spdNii_");

    system("cp ${tmpDir}/spdNii_exitcode.nii.gz ${outputRoot}ExitCode.nii.gz");
    system("cp ${tmpDir}/spdNii_dt.nii.gz ${outputRoot}DT.nii.gz");

}
