#!/usr/bin/perl -w
#
# Do some motion and optionally distortion correction with antsMotionCorr
#

my $usage = qq{
Usage: antsEddyCorrect.pl --dwi dwi.nii.gz --bvals bvals --bvecs bvecs --output-root path/to/output/root --ref-image

  Required:

  --bvals
    FSL bvals file

  --bvecs
    FSL bvecs file
  
  --dwi
    4D NIFTI image containing DWI data. Must be accompanied by bvals and bvecs. 
 
  --output-root
    Prepended onto output files

  --ref
    Reference image passed to antsMotionCorr. This should probably be either all b=0 images registered to the first, 
    or just the first b=0 volume.

  --brain-mask
    Mask used for computing displacement parameters. Not used in correction, just for diagnostics


  Optional

  --transform
    Rigid or Affine (default = Affine)

};

use strict;
use FindBin qw($Bin);
use File::Basename;
use File::Copy;
use File::Path;
use File::Spec;
use Getopt::Long;

my $dwiImage = "";
my $bvals = "";
my $bvecs = "";
my $brainMask = "";
my $outputRoot = "";
my $refImage = "";
my $transformType = "Affine";

if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

GetOptions ("dwi=s" => \$dwiImage,
	    "bvals=s" => \$bvals,
	    "bvecs=s" => \$bvecs,
	    "output-root=s" => \$outputRoot,
	    "ref=s" => \$refImage,
	    "brain-mask=s" => \$brainMask,
	    "transform=s" => \$transformType)
    or die("Error in command line arguments\n");


# Set to 1 to delete intermediate files after we're done
my $cleanup=1;

# Get the directories containing programs we need
my ($antsPath, $sysTmpDir) = @ENV{'ANTSPATH', 'TMPDIR'};

my ($outputFileRoot, $outputDir) = fileparse($outputRoot);

if ( ! -d $outputDir ) {
    mkpath($outputDir, {verbose => 0, mode => 0755}) or die "Cannot create output directory $outputDir\n\t";
}

# Directory for temporary files that is deleted later if $cleanup
my $tmpDir = "";

my $tmpDirBaseName = "${outputFileRoot}dtiproc";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = $outputDir . "/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = $sysTmpDir . "/${tmpDirBaseName}";
}

mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir\n\t";


if ( !(-f $bvals && -f $bvecs) ) {
    die("Missing bvecs or bvals\n\t");
}

# done with args

my $outputDWI = "${outputRoot}.nii.gz";
my $bvalsCorrected = "${outputRoot}.bval";
my $bvecsCorrected = "${outputRoot}.bvec";

# Bvals not changed, just copied for convenience
copy($bvals, $bvalsCorrected);

# antsMotionCorr applies same transform to all scans, not ideal because b=0 to b=0 should really 
# be rigid, but hopefully will converge on rigid solution
system("${antsPath}/antsMotionCorr -d 3 -m MI[${refImage},${dwiImage}, 1, 32, Regular, 0.25] -u 1 -t ${transformType}[0.2] -i 20 -e 1 -f 1 -s 0 -o [${outputRoot}, ${outputDWI}] --use-histogram-matching 0 --verbose");

my $mocoParams = "${outputRoot}MOCOparams.csv";

# Correct bvecs
system("${antsPath}/antsMotionCorrDiffusionDirection --bvec $bvecs --output $bvecsCorrected --moco $mocoParams --physical ${refImage}");

# Output absolute and relative displacements
system("${antsPath}/antsMotionCorrStats -x $brainMask -m $mocoParams -o ${outputRoot}AbsoluteDisplacement.csv");
system("${antsPath}/antsMotionCorrStats -x $brainMask -m $mocoParams -f 1 -o ${outputRoot}RelativeDisplacement.csv");

# cleanup
if ($cleanup) {
    
    system("rm -f ${tmpDir}/*");
    system("rmdir ${tmpDir}");
}

