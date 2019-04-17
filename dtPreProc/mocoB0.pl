#!/usr/bin/perl -w
#
# Motion-correct DWI data and bvecs for further processing
#


use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;
use File::Basename;
use Getopt::Long;

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


my $usage = qq{

  $0 
     --dwi \
     --bvecs \
     --bvals \
     --output-file \
     [options]

  
  Extracts the zero volumes from a data set, and rigidly aligns them to the first one. 

  Outputs the average b=0, aligned to the first b=0 volume. We do this (rather than using the average) because
  the first b=0 is the reference space for eddy.


  Required args

    --dwi
      4D DWI data, from which the b=0 data is extracted
  
    --bvecs
      bvecs for the data

    --bvals
      bvals for the data

    --output
      Output file, including path and ending in .nii.gz.


  Options

   --max-unweighted-b
     Specify a maximum b > 0 to include as a zero measurement. Used for schemes where the b-value of unweighted
     data is not zero. Default = 10


  Requires ANTs and Camino

};


if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

my ($dwiData, $bvecs, $bvals, $outputFile);

# Use anything <= this as a zero measurement
my $maxB = 10;


GetOptions ("dwi=s" => \$dwiData,
            "bvals=s" => \$bvals,
            "bvecs=s" => \$bvecs,
            "max-unweighted-b" => \$maxB,
	    "output-file=s" => \$outputFile
    )
    or die("Error in command line arguments\n");


my ($outputFileRoot,$outputDir) = fileparse($outputFile, (".nii.gz"));

if (! -d $outputDir ) { 
  mkpath($outputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputDir\n\t";
}

# Set to 1 to delete intermediate files after we're done
my $cleanup=1;

# Directory for temporary files that is deleted later if $cleanup
my $tmpDir = "";

my $tmpDirBaseName = "${outputFileRoot}mocoDWI";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = $outputDir . "/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = $sysTmpDir . "/${tmpDirBaseName}";
}

# Gets removed later, so check we can create this and if not, exit immediately
mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir (maybe it exists from a previous failed run)\n\t";

my $scheme = "${tmpDir}/schemefile";

# Make scheme file so that we can select b=0 data
system("fsl2scheme -bvals $bvals -bvecs $bvecs -bscale 1 -outputfile $scheme");

my $bvalString = `cat $bvals`;

chomp($bvalString);

my @bvals = split('\s+', $bvalString);

my $numB0 = 0;

foreach my $bval (@bvals) {
    if ($bval <= $maxB) {
	$numB0 = $numB0 + 1;
    }
}

if ($numB0 == 0) {
    die("No reference image for motion correction\n\t");
}

# Select b=0 data volumes
system("selectshells -inputfile $dwiData -maxbval $maxB -schemefile $scheme -outputroot ${tmpDir}/allB0");

if ($numB0 > 1) {
    system("split4dnii -inputfile ${tmpDir}/allB0.nii.gz -outputroot ${tmpDir}/b0_");
}
else {
    system("mv ${tmpDir}/allB0.nii.gz ${tmpDir}/b0_0001.nii.gz");
}

# Rigidly align all to the first one
my @b03d = `ls ${tmpDir}/b0_*.nii.gz`;

chomp(@b03d);

my $fixed = shift(@b03d);

if (scalar(@b03d) == 0) {
    print "One b=0 volume found, no registration to do\n";
    system("${antsPath}ImageMath 3 $outputFile TruncateImageIntensity ${fixed} 0 0.999 256");
    exit 0;
}

my $b0Counter = 0;

foreach my $moving (@b03d) {

    print "Registering b0 volume " . ($b0Counter + 1) . " to reference b0 \n";

    system("${antsPath}antsRegistration -d 3 -w [0,0.995] -m Mattes[ $fixed , $moving , 1, 32, Regular, 0.25] -t Rigid[0.1] -f 2x1 -s 1x0vox -c [20x20,1e-6,10] -o [ ${tmpDir}/b0ToRef_${b0Counter}, ${tmpDir}/b0ToRef_${b0Counter}_deformed.nii.gz ]");
    
    $b0Counter = $b0Counter + 1;
}

# Average
system("${antsPath}AverageImages 3 ${tmpDir}/avgB0.nii.gz 0 $fixed ${tmpDir}/b0ToRef_*_deformed.nii.gz");

# Remove outlier intensities, helps with registration later, also nicer for display
# Might need to make this an option and be a bit more aggressive in data with large signal pileup
system("${antsPath}ImageMath 3 $outputFile TruncateImageIntensity ${tmpDir}/avgB0.nii.gz 0 0.999 256");


# cleanup

if ($cleanup) {
    system("rm -f $tmpDir/*");
    system("rmdir $tmpDir");
}
