#!/usr/bin/perl -w
#
# Align DT to T1 and / or template(s)
#


use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;
use File::Basename;
use Getopt::Long;

# Get the directories containing programs we need
my ($antsPath, $sysTmpDir) = @ENV{'ANTSPATH', 'TMPDIR'};

# Defaults for options mentioned in usage
my $distCorrOutputSuffix = "NormalizedToStructural";
my $sstOutputSuffix = "NormalizedToSST";
my $groupTemplateOutputSuffix = "NormalizedToGroupTemplate";
my $standardTemplateOutputSuffix = "NormalizedToStandardTemplate";

my $backgroundMD = 7E-4;

# Save the warped tensors (can be a lot of disk space)
my $saveDeformedTensors = 0;

my $usage = qq{

  $0 
     --dt 
     --template-warp-root 
     --dist-corr-warp-root 
     --output-root 
     [options]

  Warp a DT image to template space. Assumes the existence of forward warps from DT to T1, and T1 to template. 
  Optionally, forward warps from the local template to a standard template (eg MNI152) may also be specified. 

  Warp roots should be specified such that \${root}0GenericAffine.mat exists, and optionally \${root}1Warp.nii.gz also.

  Required args

   --dt
     The dt brain image.

   --dist-corr-warp-root
     Forward warps
  
   --output-root
     Output root, including directory and file root.


  One reference image must be supplied. Additional reference images may be specified to produce DT warped to different target spaces. 
  You can specify intermediate warps without a reference image. For example:

  $0 ... --dist-corr-warp-root /path/to/distcorr --template-warp-root /path/to/antsctWarp \
         --template localTemplate.nii.gz 

  This uses the distortion correction to create a composite warp to template space, but only produces output in the template space. Alternatively,

  $0 ... --dist-corr-warp-root /path/to/distcorr --template-warp-root /path/to/antsctWarp \
         --structural-image t1.nii.gz \
         --template localTemplate.nii.gz 

  This produces output in both the T1 and the template space.

  Options  

   --mask
     Brain mask for the tensor in the native space. The tensor is then masked with this image and background voxels are
     filled with isotropic tensors with a constant MD (default 7E-4). This helps avoid interpolation artifacts where tensor 
     components are interpolated with zero in the log space. Using the mask helps avoid the ring of high diffusivity around 
     the edge of the brain. However, it means that background MD / RD / AD values will be non-zero, unless masks are also 
     supplied for the template images.

   --background-md 
     Mean diffusion for background voxels. This probably doesn't need changing unless your diffusion tensors do not 
     have the usual units of mm^2 / s. Only used if the mask is specified (default = $backgroundMD).

   --save-tensors
     Save the deformed tensor(s) in the target space(s) (default = ${saveDeformedTensors}).

   --structural
     The structural image (usually T1) that is the reference image for the intra-subject distortion correction. 
     Required to resample the image into the structural image space.

   --structural-mask
     Brain mask in the structural space. Used to set the background tensors to zero after warping.

   --structural-output-suffix
     String appended to the output root to denote images in the structural space (default = $distCorrOutputSuffix).

   --single-subject-template
     The SST reference image for this subject. Required to ressample the DT into single-subject template space.

   --single-subject-template-mask
     A brain mask in the SST space. Used to set the background tensors to zero after warping.

   --single-subject-template-warp-root
     Forward warps defining a transform from the structural space to the SST.

   --single-subject-template-output-suffix
     String appended to the output root to denote images in the SST space, where the warp is a composition of
     diffusion -> structural -> template (default = $sstOutputSuffix). Should not contain a file extension.

   --group-template-warp-root 
     Forward warps defining a transform from the SST to the group template. 

   --group-template
     Group space reference image. Required to resample the DT into the group template space.

   --group-template-mask
     A brain mask in the group template space. Used to set the background tensors to zero after warping.

   --group-template-output-suffix
     String appended to the output root to denote images in the group template space, where the warp is a composition of
     diffusion -> structural -> SST -> group template (default = $groupTemplateOutputSuffix). Should not contain a file extension.

   --standard-template-warp-root 
     Forward warps defining a transform from the group template to the standard space. If this is specified, the standard template 
     image must also be specified.

   --standard-template
     Standard space reference image. Required to resample the DT into the standard template space.

   --standard-template-mask
     A brain mask in the standard template space. Used to set the background tensors to zero after warping.

   --standard-template-output-suffix
     String appended to the output root to denote images in the standard template space, where the warp is a composition of
     diffusion -> structural -> SST -> group template -> standard template (default = $standardTemplateOutputSuffix). 
     Should not contain a file extension.


  Output DT + scalar metrics are computed for each reference space provided:

    DT - diffusion tensors
    FA - fractional anisotropy
    RD - radial diffusivity (L2 + L3) / 2
    AD - Axial diffusivity (L1)
    MD - Mean diffusivity (L1 + L2 + L3) / 3 = (AD + 2*RD) / 3
    RGB - Tensor principal direction color image modulated by FA
  

  Requires ANTs

};


if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

# Require all of these
my ($dt, $outputRoot, $distCorrWarpRoot);

# Options have default settings
my $dtMask = "";
my $structuralImage = "";
my $structuralMask = ""; 

my $sst = "";
my $sstWarpRoot = "";
my $sstMask = "";

my $groupTemplate = "";
my $groupTemplateWarpRoot = "";
my $groupTemplateMask = "";


my $standardTemplate = "";
my $standardTemplateWarpRoot = "";
my $standardTemplateMask = "";


GetOptions ("dt=s" => \$dt,
	    "mask=s" => \$dtMask,
	    "dist-corr-warp-root=s" => \$distCorrWarpRoot,
	    "structural=s" => \$structuralImage,
	    "structural-mask=s" => \$structuralMask,
	    "single-subject-template=s" => \$sst,
	    "single-subject-template-mask=s" => \$sstMask,
	    "single-subject-template-warp-root=s" => \$sstWarpRoot,
	    "group-template=s" => \$groupTemplate,
	    "group-template-warp-root=s" =>\$groupTemplateWarpRoot,
	    "group-template-mask=s" => \$groupTemplateMask,
	    "standard-template=s" => \$standardTemplate,
	    "standard-template-warp-root=s" =>\$standardTemplateWarpRoot,
	    "standard-template-mask=s" => \$standardTemplateMask,
	    "background-md=f" => \$backgroundMD,
	    "single-subject-template-output-suffix=s" => \$sstOutputSuffix,
	    "group-template-output-suffix=s" => \$groupTemplateOutputSuffix,
	    "standard-template-output-suffix=s" => \$standardTemplateOutputSuffix,
	    "output-root=s" => \$outputRoot,
	    "save-tensors=i" => \$saveDeformedTensors

    )
    or die("Error in command line arguments\n");


my ($outputFileRoot,$outputDir) = fileparse($outputRoot);

# Check input
if (! -f $dt) {
    die("Tensor input $dt does not exist");
}
if (! (-f $structuralImage || -f $sst || -f $groupTemplate || -f $standardTemplate)) {
    die("At least one reference image is required");
}

if (! -d $outputDir ) { 
  mkpath($outputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputDir\n\t";
}

# Set to 1 to delete intermediate files after we're done
# Has no long term effect if using qsub since files get cleaned up anyhow
my $cleanup=1;

# Directory for temporary files that is deleted later if $cleanup
my $tmpDir = "";

my $tmpDirBaseName = "${outputFileRoot}dtToTemplate";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = $outputDir . "/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = $sysTmpDir . "/${tmpDirBaseName}";
}

# Gets removed later, so check we can create this and if not, exit immediately
mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir (maybe it exists from a previous failed run)\n\t";


# Construct warp strings for antsApplyTransforms
my $distCorrWarpString = getWarpString($distCorrWarpRoot);
my $sstWarpString = "";
my $groupTemplateWarpString = "";
my $standardTemplateWarpString = "";


# Get all warps needed to take us to a particular target space
# getWarpString will die if the warps don't exist
if (-f $standardTemplate) {
    $standardTemplateWarpString = getWarpString($standardTemplateWarpRoot);
    $groupTemplateWarpString = getWarpString($groupTemplateWarpRoot);
    $sstWarpString = getWarpString($sstWarpRoot);
}

if (-f $groupTemplate) {
    $groupTemplateWarpString = getWarpString($groupTemplateWarpRoot);
    $sstWarpString = getWarpString($sstWarpRoot);
}

if (-f $sst) {
    $sstWarpString = getWarpString($sstWarpRoot);
}

# Create output in the target spaces

# File naming convention: ${outputRoot}DT${suffix}.nii.gz

# Scalar images will be the same thing but "DT" is replaced with "FA", etc

my $dtMoving = $dt;

if (-f $dtMask) {
    $dtMoving = "${tmpDir}/dtMasked.nii.gz";
    system("${antsPath}ImageMath 3 $dtMoving TensorMask $dt $dtMask $backgroundMD")
}

if (-f $structuralImage) {
    print "\nCreating output in space of $structuralImage\n";

    my $dtDistCorr = "${tmpDir}/dt${distCorrOutputSuffix}.nii.gz";

    warpAndReorientDT($dtMoving, $structuralImage, $distCorrWarpString, $structuralMask, $backgroundMD, $dtDistCorr);

    createScalarImages($dtDistCorr, $outputRoot, $distCorrOutputSuffix);

    if ($saveDeformedTensors) {
	system("mv $dtDistCorr ${outputRoot}DT${distCorrOutputSuffix}.nii.gz");
    }

}

if (-f $sst) {
    print "\nCreating output in space of $sst\n";

    my $dtSST = "${tmpDir}/dt${sstOutputSuffix}.nii.gz";

    warpAndReorientDT($dtMoving, $sst, "$sstWarpString $distCorrWarpString", $sstMask, $backgroundMD, $dtSST);

    createScalarImages($dtSST, $outputRoot, $sstOutputSuffix);

    if ($saveDeformedTensors) {
	system("mv $dtSST ${outputRoot}DT${sstOutputSuffix}.nii.gz");
    }

}

if (-f $groupTemplate) {
    print "\nCreating output in space of $groupTemplate\n";

    my $dtGroupTemplate = "${tmpDir}/dt${groupTemplateOutputSuffix}.nii.gz";

    warpAndReorientDT($dtMoving, $groupTemplate, "$groupTemplateWarpString $sstWarpString $distCorrWarpString", 
		      $groupTemplateMask, $backgroundMD, $dtGroupTemplate);

    createScalarImages($dtGroupTemplate, $outputRoot, $groupTemplateOutputSuffix);

    if ($saveDeformedTensors) {
	system("mv $dtGroupTemplate ${outputRoot}DT${groupTemplateOutputSuffix}.nii.gz");
    }

}

if (-f $standardTemplate) {
    print "\nCreating output in space of $standardTemplate\n";

    my $dtStandardTemplate = "${tmpDir}/dt${standardTemplateOutputSuffix}.nii.gz";

    warpAndReorientDT($dtMoving, $standardTemplate, 
		      "$standardTemplateWarpString $groupTemplateWarpString $sstWarpString $distCorrWarpString", 
		      $standardTemplateMask, $backgroundMD, $dtStandardTemplate);

    createScalarImages($dtStandardTemplate, $outputRoot, $standardTemplateOutputSuffix);

    if ($saveDeformedTensors) {
	system("mv $dtStandardTemplate ${outputRoot}DT${standardTemplateOutputSuffix}.nii.gz");
    }

}


# cleanup

if ($cleanup) {
    system("rm -f ${tmpDir}/*");
    system("rmdir $tmpDir");
}


#
# warpAndReorientDT($dt, $refImage, $warpString, $refMask, $backgroundMD, $outputFile);
#
# Requires $tmpDir to be defined in scope
# 
sub warpAndReorientDT {

    # where $warpString is something like "-t t1ToTemplate1Warp.nii.gz -t t1ToTemplate0GenericAffine.mat -t dtToTemplate1Warp.nii.gz -t dtToTemplate0GenericAffine.mat"
    my ($inputDT, $ref, $warpString, $refMask, $backgroundMD, $outputFile) = @_;

    my $dtTmp = "${tmpDir}/dtWarpedButNotReoriented.nii.gz";

    my $dtDF = "${tmpDir}/dtCombinedWarpField.nii.gz";
    
    system("${antsPath}antsApplyTransforms -d 3 -i $inputDT -e 2 -r $ref $warpString -o $dtTmp --verbose --float");
    system("${antsPath}antsApplyTransforms -d 3 -i $inputDT -e 2 -r $ref $warpString -o [${dtDF}, 1] --verbose --float");

    if (-f $refMask) {
	system("${antsPath}ImageMath 3 $dtTmp TensorMask $dtTmp $refMask 0");
    }

    system("${antsPath}ReorientTensorImage 3 $dtTmp $outputFile $dtDF");

    # Clean up temp files so that they can't get accidentally reused if this method is call again
    system("rm $dtTmp $dtDF");
    
}

#
# $warps = getWarpString($root)
#
# Looks for ${root}[1Warp.nii.gz | 0GenericAffine.mat]
#
sub getWarpString {

    my ($warpRoot) = @_;
    
    my $warpString = "-t ${warpRoot}0GenericAffine.mat";
    
    
    if (! -f "${warpRoot}0GenericAffine.mat") {
	die "\nNo affine transform matching ${warpRoot}0GenericAffine.mat\n";
    }

    if (-f "${warpRoot}1Warp.nii.gz") {
	$warpString = "-t ${warpRoot}1Warp.nii.gz $warpString";
    }
  
    return $warpString;

}

#
# createScalarImages($dt, $outputRoot, $outputSuffix)
#
# Computes ${outputRoot}[FA, MD, RD, AD, RGB]${outputSuffix}
#
sub createScalarImages {

    my ($dt, $outputRoot, $outputSuffix) = @_;

    system("${antsPath}ImageMath 3 ${outputRoot}FA${outputSuffix}.nii.gz TensorFA $dt");
    system("${antsPath}ImageMath 3 ${outputRoot}MD${outputSuffix}.nii.gz TensorMeanDiffusion $dt");
    system("${antsPath}ImageMath 3 ${outputRoot}RD${outputSuffix}.nii.gz TensorRadialDiffusion $dt");
    system("${antsPath}ImageMath 3 ${outputRoot}AD${outputSuffix}.nii.gz TensorAxialDiffusion $dt");
    system("${antsPath}ImageMath 3 ${outputRoot}RGB${outputSuffix}.nii.gz TensorColor $dt");
    
}
