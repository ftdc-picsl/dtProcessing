#!/usr/bin/perl -w
#
# Average deformed DTs in SST space
#


use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;
use File::Basename;
use Getopt::Long;

my $antsPath = "/data/grossman/pipedream2018/bin/ants/bin/";

# Get the directories containing programs we need
my ($sysTmpDir) = $ENV{'TMPDIR'};

# Defaults for options mentioned in usage

my $backgroundMD = 7E-4;


my $usage = qq{

  $0 
     --dt-norm-dir
     --subject
     [options]

  Only works if you have saved the DT in SST space. Otherwise use sstAverageScalars.pl

  Required args

   --subject
     Subject ID.

   --dt-norm-dir
     Path where the DTs to be averaged are in dtNormDir/subject/tp/subject_tp_DTNormalizedToSST.nii.gz
     Output is in dtNormDir/subject/singleSubjectTemplate.

  Options  

   --mask
     Brain mask for the SST.

   --save-tensor
     If 1, save the diffusion tensor (default = 0).

  Scalar metrics are computed after averaging

    DT - diffusion tensors
    FA - fractional anisotropy
    MD - Mean diffusivity (L1 + L2 + L3) / 3 = (AD + 2*RD) / 3
    RGB - Tensor principal direction color image modulated by FA


};


if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

# Require all of these
my ($subject, $dtNormDir);

# Options have default settings
my $sstMask = "";
my $saveTensor = 0;

GetOptions ("subject=s" => \$subject,
	    "dt-norm-dir=s" => \$dtNormDir,
	    "mask=s" => \$sstMask,
	    "save-tensor=i" => \$saveTensor
    )
    or die("Error in command line arguments\n");


# Set to 1 to delete intermediate files after we're done
# Has no long term effect if using qsub since files get cleaned up anyhow
my $cleanup=1;


my $outputDir = "${dtNormDir}/${subject}/singleSubjectTemplate";

if (! -d $outputDir ) { 
  mkpath($outputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputDir\n\t";
}


my $tmpDir = "";

my $tmpDirBaseName = "${subject}_sstAverageDT";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = $outputDir . "/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = $sysTmpDir . "/${tmpDirBaseName}";
}

# Gets removed later, so check we can create this and if not, exit immediately
mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir (maybe it exists from a previous failed run)\n\t";


# Check input

my @tps = `ls ${dtNormDir}/$subject | grep -v singleSubjectTemplate`;

chomp @tps;

my $imageString = "";

foreach my $tp (@tps) {
    $imageString = $imageString . " ${dtNormDir}/${subject}/${tp}/dtNorm/${subject}_${tp}_DTNormalizedToSST.nii.gz";
}


my $avgDT = "${tmpDir}/${subject}_AverageDT.nii.gz";

my $cmd = "${antsPath}AverageTensorImages 3 $avgDT 0 $imageString";

system("$cmd");

# Check average is successfully computed
if (! -f $avgDT) {
    die("\n  Error computing average DT from $imageString");
}


my $outputRoot = "${outputDir}/${subject}_Average";

createScalarImages($avgDT, $outputRoot, "");

if ($saveTensor) {
    system("cp $avgDT ${outputRoot}DT.nii.gz");
}

# cleanup

if ($cleanup) {
    system("rm -f ${tmpDir}/*");
    system("rmdir $tmpDir");
}



#
# createScalarImages($dt, $outputRoot, $outputSuffix)
#
# Computes ${outputRoot}[FA, MD, RGB]${outputSuffix}
#
sub createScalarImages {

    my ($dt, $outputRoot, $outputSuffix) = @_;

    system("${antsPath}ImageMath 3 ${outputRoot}FA${outputSuffix}.nii.gz TensorFA $dt");
    system("${antsPath}ImageMath 3 ${outputRoot}MD${outputSuffix}.nii.gz TensorMeanDiffusion $dt");
    system("${antsPath}ImageMath 3 ${outputRoot}RGB${outputSuffix}.nii.gz TensorColor $dt");
    
}
