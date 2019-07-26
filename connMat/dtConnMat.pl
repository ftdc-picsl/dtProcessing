#!/usr/bin/perl -w
#
# Do tractography in the DT, then evaluate connectivity in structural space
#


use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;
use File::Basename;
use Getopt::Long;

# Get the directories containing programs we need
my ($antsPath, $sysTmpDir) = @ENV{'ANTSPATH', 'TMPDIR'};


my $usage = qq{

  $0 
     --dt 
     --mask
     --reference-image
     --dist-corr-warp-root | --dt-to-reference-tract-warp 
     --label-image
     --label-def
     --exclusion-image
     --output-root 
     [options]

  Does tractography in the DT space and processes results in another space.

  Warp root should be specified such that \${root}0GenericAffine.mat exists, and optionally \${root}1Warp.nii.gz also.
  More complex warps (eg to a longitudinal group template) should pass a single composed warp field.

  The FA image is upsampled to 1mm isotropic and thresholded to make a seed mask.


  Required args

   --dt
     The dt brain image.

   --mask 
     The dt brain mask.

   --reference-image     
     The reference image, eg the subject's T1 image.
 
   --dist-corr-warp-root 
     Root of warps, where the reference image is the fixed image. This is used when there is a single Warp / Affine pair 
     of transforms, usually a distortion correction mapping the the DWI to a T1 image. For more complex transforms, use
      --dt-to-reference-tract-warp to pass a single warp field to use to map streamlines. Note that the order of transforms for 
     warping point sets is reversed with respect to the order required for warping images.

   --dt-to-reference-tract-warp
     A composed warp for transforming point sets into the reference space. 

   --label-image
     Image containing target labels in the reference space, these are the nodes of the connectivity graph.

   --label-def
     CSV file containing the label IDs and names of the labels to be included in the graph. Labels in the image that are 
     not in this list are ignored in the connectivity computation.

   --exclusion-image
     A mask or probability image of voxels where tracts should be truncated in the reference space. This would usually be a CSF mask. 
     A probability or label image may be passed here, the image will be thresholded at 0.5. This truncates tracts in the reference space
     and does not affect the generation of streamlines. The tractography algorithm is constrained to voxels included in the DT mask passed 
     with --mask. 

   --output-root
     Output root, including directory and file root.

  Options  

   --count-longest-path
     Compute the longest connection between end points, rather than the shortest. If you do this, you need 
     to constrain tracking somehow (eg, with a WM mask) or you may get bad results because of tracking
     continuing through GM into unconnected regions.

   --seed-spacing 
     Spacing in mm of seed points. The DT image is resampled to be isotropic at this resolution, FA is then 
     computed and the thresholds are applied to generate a seed image (default = 1).

   --seed-fa-thresh
     Threshold at which to define seed points (default = 0.25).

   --curve-thresh
      Maximum curvature, in degrees, over which streamlines are allowed to bend over 5mm intervals (default = 90).

   --compute-scalars
     Compute the median tract (FA, RD, AD, MD) in addition to the streamline counts (default = 0).

   --dt-to-reference-image-warp
     A composed warp for transforming images sets into the reference space. This is used to warp the DT into the reference 
     space, after which FA, MD, etc are sampled. Required if using if using --dt-to-reference-tract-warp and scalars are
     to be computed.


  Requires ANTs, ANTsR and Camino

};


if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

if (!$antsPath || ! -f "${antsPath}antsRegistration") {
    die("Script requires ANTSPATH to be defined");
}

my $track = `which track`;

chomp $track;

if (! -f $track) {
    die("Script requires Camino");
}

# Require all of these
my ($dt, $fa, $mask, $outputRoot, $labelImage, $labelDef, $exclusion, $referenceImage);

# Require one of these
my $distCorrWarpRoot = "";

my $composedTractWarp = "";
my $composedImageWarp = "";

# Options have default settings
my $seedMinFA = 0.25;
my $curveThreshDeg = "90";
my $seedSpacing = 1;
my $computeScalars = 0;
my $countLongestPath = 0;

GetOptions ("dt=s" => \$dt,
	    "exclusion-image=s" => \$exclusion,
	    "dist-corr-warp-root=s" => \$distCorrWarpRoot,
	    "dt-to-reference-tract-warp=s" => \$composedTractWarp,
	    "dt-to-reference-image-warp=s" => \$composedImageWarp,
	    "label-image=s" => \$labelImage,
	    "label-def=s" => \$labelDef,
	    "mask=s" => \$mask,
	    "output-root=s" => \$outputRoot,
	    "reference-image=s" => \$referenceImage,
	    "seed-spacing=f" => \$seedSpacing,
	    "seed-fa-thresh=f" => \$seedMinFA,
	    "curve-thresh=f" => \$curveThreshDeg,
	    "compute-scalars=i" => \$computeScalars,
	    "count-longest-path=i" => \$countLongestPath
    )
    or die("Error in command line arguments\n");


if (-f $composedTractWarp) {
    if ($computeScalars) {
	if (! -f $composedImageWarp) {
	    die("\n  Missing required forward warp to transform DT to reference space for scalar matrix computation (see usage for --dt-to-reference-image-warp)\n");
	}
    }
}
else {
    if (! (-f "${distCorrWarpRoot}1InverseWarp.nii.gz" && -f "${distCorrWarpRoot}0GenericAffine.mat") ) {
	die("\n  Missing required transform to reference space\n");
    }
}

# Some other settings we hard code

# Discard anything shorter than this (mm)
my $minTractLength = 10;

# For historical reasons, Camino doesn't support an explicit brain mask with DT tractography. Tracking is restricted to the brain by using 
# the brain mask as an anisotropy image. To incorporate an FA threshold on the tracking, you would need to threshold FA and multiply it
# by the brain mask.

my ($outputFileRoot, $outputDir) = fileparse($outputRoot);


if (! -d $outputDir ) { 
  mkpath($outputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputDir\n\t";
}

# Set to 1 to delete intermediate files after we're done
my $cleanup=0;

# Directory for temporary files that is deleted later if $cleanup
my $tmpDir = "";

my $tmpDirBaseName = "${outputFileRoot}dtConnMat";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = $outputDir . "/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = $sysTmpDir . "/${tmpDirBaseName}";
}

# Gets removed later, so check we can create this and if not, exit immediately
mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir (maybe it exists from a previous failed run)\n\t";

print "\nCreating seed image at $seedSpacing mm isotropic resolution, min FA $seedMinFA\n";

# Resample the DT and get FA
my $seedSpaceRefImage = "${tmpDir}/faResample.nii.gz";

system("${antsPath}ImageMath 3 ${tmpDir}/fa.nii.gz TensorFA $dt");

# Resampled FA just used to provide a reference space
system("${antsPath}ResampleImageBySpacing 3 ${tmpDir}/fa.nii.gz $seedSpaceRefImage $seedSpacing $seedSpacing $seedSpacing 0 0 0");

# Now resample the DT, and compute FA from that
system("${antsPath}antsApplyTransforms -d 3 -e 2 -i $dt -o ${tmpDir}/dtResample.nii.gz -r $seedSpaceRefImage");

system("${antsPath}ImageMath 3 ${tmpDir}/faForSeeds.nii.gz TensorFA ${tmpDir}/dtResample.nii.gz");

my $seedImage = "${tmpDir}/seedMask.nii.gz";

system("${antsPath}ThresholdImage 3 ${tmpDir}/faForSeeds.nii.gz $seedImage $seedMinFA Inf");

print "\nTracking and warping streamlines to reference space \n";

# Track, and output lines in structural space

# Get appropriate point set warp
my $warpString = "";

if (-f $composedTractWarp) {
    $warpString = "--composed-inverse-warp $composedTractWarp";
}
else {
    $warpString = "--dist-corr-inverse-warp=${distCorrWarpRoot}1InverseWarp.nii.gz --dist-corr-affine=${distCorrWarpRoot}0GenericAffine.mat";
}

system("track -inputmodel dt -inputfile $dt -anisfile $mask -anisthresh 0.5 -curvethresh $curveThreshDeg -seedfile $seedImage -tracker euler -interpolator linear -stepsize 0.5 -silent | ${Bin}/warpTracts.R --input-file=stdin --output-root=${tmpDir}/ --seed-image=$seedImage --reference-image=$referenceImage $warpString");

my $allTracts = "${tmpDir}/TractsDeformed.Bfloat";

# Threshold exclusion image, in case it is a probability
system("${antsPath}ThresholdImage 3 $exclusion ${tmpDir}/exclusion.nii.gz 0.5 Inf");

# For debug purposes, compute ACM before and after exclusion etc

# Produce a global ACM of all tracts 
system("procstreamlines -inputfile $allTracts -outputacm -header $referenceImage -outputroot ${tmpDir}/allTracts_ -silent");

# Produce a global ACM of all tracts with some junk taken out 
system("procstreamlines -inputfile $allTracts -mintractlength $minTractLength -exclusionfile ${tmpDir}/exclusion.nii.gz -truncateinexclusion -outputacm -header $referenceImage -outputroot ${tmpDir}/allTractsWithExclusion_ -silent");

# Filter and build connectivity matrix, plus other useful diagnostics
my $endpointTracts = "${tmpDir}/endpointTracts.Bfloat";
my $graphTracts = "${tmpDir}/graphTracts.Bfloat";

print "\nComputing connectivity matrix\n";

my $cmd = "procstreamlines -inputfile $allTracts -mintractlength $minTractLength -exclusionfile ${tmpDir}/exclusion.nii.gz -truncateinexclusion -endpointfile $labelImage -header $referenceImage -outputfile $endpointTracts -silent";

if ($countLongestPath) {
    $cmd = $cmd . " -allowmultitargets";
}

system($cmd);

# Because raw streamline file is potentially very large, remove as soon as we're done with it
system("rm -f ${allTracts}");

# Further remove tracts below min tract length after exclusion / end points
system("procstreamlines -inputfile $endpointTracts -mintractlength $minTractLength -header $referenceImage -outputfile $graphTracts -silent");

system("rm -f ${endpointTracts}");

# Produce a global ACM of all tracts that contribute to the connectivity matrix
system("procstreamlines -inputfile $graphTracts -outputacm -header $referenceImage -outputroot ${tmpDir}/graphTracts_ -silent");

$cmd = "conmat -inputfile $graphTracts -outputroot ${tmpDir}/conmat_ -targetfile $labelImage -targetnamefile $labelDef -tractstat length";

if ($countLongestPath) {
    $cmd = $cmd . " -countlongestpath";
}

system($cmd);

# Move results to output location
system("cp ${tmpDir}/conmat_sc.csv ${outputRoot}sc.csv");
system("cp ${tmpDir}/conmat_ts.csv ${outputRoot}MeanTractLength.csv");
system("cp ${tmpDir}/graphTracts_acm_sc.nii.gz ${outputRoot}GraphTractsACM.nii.gz");
system("cp ${tmpDir}/allTracts_acm_sc.nii.gz ${outputRoot}AllTractsACM.nii.gz");
system("cp ${tmpDir}/allTractsWithExclusion_acm_sc.nii.gz ${outputRoot}AllTractsACMWithExclusion.nii.gz");
system("cp ${tmpDir}/SeedsDeformed.nii.gz ${outputRoot}SeedDensityDeformed.nii.gz");

if ($computeScalars) {

    # Resample DT to Structural space
    print "\nComputing diffusion scalar connectivity matrices \n";

    my $dtWarpString = "";

    if (-f $composedImageWarp) {
	$dtWarpString = "-t $composedImageWarp";
    }
    else {
	$dtWarpString = "-t ${distCorrWarpRoot}1Warp.nii.gz -t ${distCorrWarpRoot}0GenericAffine.mat";
    }

    system("${antsPath}antsApplyTransforms -d 3 -e 2 -i $dt -o ${tmpDir}/dtStructural.nii.gz $dtWarpString -r $labelImage ");

    system("${antsPath}ImageMath 3 ${tmpDir}/dtStructuralFA.nii.gz TensorFA ${tmpDir}/dtStructural.nii.gz");
    system("${antsPath}ImageMath 3 ${tmpDir}/dtStructuralMD.nii.gz TensorMeanDiffusion ${tmpDir}/dtStructural.nii.gz");
    system("${antsPath}ImageMath 3 ${tmpDir}/dtStructuralAD.nii.gz TensorAxialDiffusion ${tmpDir}/dtStructural.nii.gz"); 
    system("${antsPath}ImageMath 3 ${tmpDir}/dtStructuralRD.nii.gz TensorRadialDiffusion ${tmpDir}/dtStructural.nii.gz");

    my @scalars = qw(FA MD AD RD);

    foreach my $scalar (@scalars) {
	$cmd = "conmat -inputfile $graphTracts -outputroot ${tmpDir}/conmat${scalar}_ -targetfile $labelImage -targetnamefile $labelDef -tractstat median -scalarfile ${tmpDir}/dtStructural${scalar}.nii.gz";

	if ($countLongestPath) {
	    $cmd = $cmd . " -countlongestpath";
	}
	
	system($cmd);
	
	system("cp ${tmpDir}/conmat${scalar}_ts.csv ${outputRoot}MeanTractMedian${scalar}.csv");
    }

}


# cleanup

if ($cleanup) {
    system("rm -f ${tmpDir}/*");
    system("rmdir $tmpDir");
}
