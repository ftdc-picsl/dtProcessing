#!/usr/bin/perl -w
#
# Align DWI to intra-session structural image
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


my $usage = qq{

  $0 
     --structural 
     --mask
     --dwi b0.nii.gz 
     --output-root

  Corrects a DWI image (usually b0) to the structural image (usually T1). It is assumed, but not required, that the DWI image is not skull stripped.

  Required args

   --structural
     The structural brain image. Should be brain extracted.

   --dwi
     The DWI image, usually b0. Does not have to be brain extracted, we use a registration mask + assumption of small deformation to avoid problems.
  
   --mask
     A brain mask in the structural space. This will be transferred to the DWI space.

   --output-root
     Output root, including directory and file root.


  Output prefixed with the output root, then:

   0GenericAffine.mat 
   1Warp.nii.gz
   1InverseWarp.nii.gz  - distortion correction, forward warps deform DWI to T1 space

   DWIDeformed.nii.gz   - DWI image, deformed to T1 space

   DWIBrainMask.nii.gz  - T1 brain mask warped to DWI space  

  
  Requires ANTs

};


if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

my ($structuralImage, $brainMask, $dwi, $outputRoot);


GetOptions ("dwi=s" => \$dwi,
	    "structural=s" => \$structuralImage,  
            "mask=s" => \$brainMask,
	    "output-root=s" => \$outputRoot
    )
    or die("Error in command line arguments\n");


my ($outputFileRoot,$outputDir) = fileparse($outputRoot);

if (! -d $outputDir ) { 
  mkpath($outputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputDir\n\t";
}

# Set to 1 to delete intermediate files after we're done
# Has no effect if using qsub since files get cleaned up anyhow
my $cleanup=1;

# Directory for temporary files that is deleted later if $cleanup
my $tmpDir = "";

my $tmpDirBaseName = "${outputFileRoot}dwiDistCorr";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = $outputDir . "/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = $sysTmpDir . "/${tmpDirBaseName}";
}

# Gets removed later, so check we can create this and if not, exit immediately
mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir (maybe it exists from a previous failed run)\n\t";


my $distCorrRegMask="${tmpDir}/distcorrRegistrationMask.nii.gz";

system("${antsPath}ImageMath 3 ${distCorrRegMask} MD $brainMask 4");

#
# Could apply some restrict deformation here eg -g 0.25x1x0.25
# Have to be careful because this vector is applied in the voxel space, so to do this you have to know the phase encoding direction. Additionally, we must
# assume rigid rotation is small, and allow some off axis deformation (needed even without motion to prevent numerical problems in ANTs, eg 0.1x1x0.1)
#
# BSplineSyN[ gradientStep, meshSpacingForUpdate, meshSpacingForTotal, splineOrder]
#
# The mesh spacing should be a single number, eg 40 for 40mm. 
# Alternatively, you can set size as a vector, eg 10x10x10 means calculate a knot spacing such that there are 10 control points in each direction.
# The first number is for the update field, the second for the total field.
# Set the second number to nonzero for some regularization.
#
# Another idea is to use different mesh sizes for each dimension, eg 4x6x6. This is tricky because specifying mesh size means the spacing depends on the
# dimensions of the input image. But the goal is to accomplish the same thing as before, make the deformation field more detailed in the phase-encode direction.
#

#
# If there are problems with aligning the edge of the brain, maybe try a tighter brain mask for the deformable stage.
#

my $fixed = $structuralImage;
my $moving = $dwi;

system("${antsPath}antsRegistration -d 3 -u 0 -o [ ${outputRoot}, ${outputRoot}DWIDeformed.nii.gz ] -t Rigid[0.1] -m Mattes[ $fixed, $moving, 1, 16, Regular, 0.25] -c [30x30x0, 1e-7,10] -f 3x2x1 -s 2x1x0vox -x $distCorrRegMask -t BSplineSyN[0.2, 40, 0, 3] -m CC[ $fixed, $moving, 1, 4] -c [20x20x0,1e-7,10] -f 2x2x1 -s 1x0x0vox -v 1 --float");

system("${antsPath}antsApplyTransforms -d 3 -i $brainMask -r $dwi -t [ ${outputRoot}0GenericAffine.mat , 1] -t ${outputRoot}1InverseWarp.nii.gz -n NearestNeighbor -o ${outputRoot}DWIBrainMask.nii.gz --float"); 



# cleanup

if ($cleanup) {
    system("rm -rf $tmpDir");
}
