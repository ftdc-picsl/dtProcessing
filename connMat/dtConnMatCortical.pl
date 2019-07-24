#!/usr/bin/perl -w
#
# Wrapper script for calling dtConnMat.pl for some or all of a subject's data
#

use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;
use File::Basename;
use Getopt::Long;

my ($antsPath, $sysTmpDir) = @ENV{'ANTSPATH', 'TMPDIR'};


# CSV file containing JLF cortical label definitions
my $jlfCorticalLabelDef = "${Bin}/mindBoggleCorticalGraphNodes.csv";

# CSV file containing JLF WM label definitions
my $jlfWMLabelDef = "${Bin}/mindBoggleWMLabels.csv";

# These are what are used to label the cortical mask from JLF, and become the graph nodes
# Default to the JLF labels themselves
my $graphCorticalLabelSystem = "PG_antsLabelFusionLabels";
my $graphCorticalLabelDef = $jlfCorticalLabelDef;

my $usage = qq{

  $0  
      --subject
      --timepoint
      --antsct-base-dir | --antslongct-base-dir
      --dt-base-dir
      --output-base-dir
      [ options ]

  Required args:

   --subject
     Subject ID

   --timepoint  
     Timepoint to process.

   --antsct-base-dir | --antslongct-base-dir  
     Base ants[long]CT dir for T1 data. If longitudinal, set --longitudinal-target.

   --dt-base-dir
     Base DTI dir. There should be DT data for the time point(s) to be processed. For longitudinal data,
     there should also be FA images in the single-subject template space.

   --output-base-dir
     Base output directory. Output is organized under this as subj/tp/connMat_corticalLabelSystem.


  Options:

   --longitudinal-target
     Required only for longitudinal data. Either "session" to evaluate connectivity in the intra-session 
     T1, or "sst", to do it in the SST space. If "session", graph nodes must be defined for each time point.
  
   --cortical-label-system
     System (without .nii.gz) for an image in the timepout antsCT output directory containing cortical labels. 
     By default, this is the JLF labels. If this option is specified, the list of cortical label IDs and names 
     must be provided with the --cortical-label-def option. The image may contain non-cortical labels, but only 
     the cortical labels will be used for the graph (default = $graphCorticalLabelSystem).

   --cortical-label-def
     A CSV file containing "Label.ID,Label.Name", with cortical labels only (default = $graphCorticalLabelSystem).

  Some pre-processing is done at run time so it is best to run this script from a qlogin session, or qsub it. 
  You may need to set CAMINO_HEAP_SIZE first to allocate enough RAM for the creation of the nodes.

  The target T1 image / SST should be labeled using the pseudo-geodesic JLF. The JLF is always used to define
  the white matter and cortical masks, after which the cortical labels may be replaced with a custom label set.
  This allows multiple parcellations to be used with a constant WM mask.

  DT data is read from

    \${dt-base-dir}/subj/tp/

  The DTI preprocessing pipeline should have been run, such that dt/ and distCorr/ exist; 
  these will be used to do the tracking and transfer the results to T1 space.

  
  Output:

    Prefixed with outputBaseDir/subj/tp/connMat_corticalLabelSystem/subj_tp_ :

      GraphNodes.nii.gz     
 
    Cortical nodes, masked to exclude high FA without introducing holes.


      ExclusionMask.nii.gz

    A mask defined by the dilation of the JLF WM segmentation and adjoining regions of high FA.
    This constrains tracking to white matter, but still allows lines to reach the nodes.


      CorticalMask.nii.gz

    A mask of all voxels containing the nodes.


      CorticalMaskEdits.nii.gz

    Edits made to the cortical mask from JLF. This means voxels newly defined as cortex.


      LabeledTractInclusionMask.nii.gz

    A label image containing all voxels that tracts can traverse without being terminated, including
    white matter and the cortical labels.


      CorticalLabelDiffMask.nii.gz
    
    A mask that is 1 wherever the cortical label in the graph nodes differs from the input cortical label 
    image. This can be used to see differences introduced by propagating labels into the JLF cortical mask.
    
};

if ($#ARGV < 0) {
    print $usage;
    exit 1;
}


# Input base dir for DT data
my $dtBaseDir = "";

# For T1 brain and other useful stuff
my $antsCTBaseDir = "";

# output to $outputBaseDir/subj/tp/connMat
my $outputBaseDir = "";

my $subject = "";
my $timepoint = "";

# If empty string, expect cross-sectional antsct input
my $longitudinalTarget = "";

GetOptions ("subject=s" => \$subject,
	    "timepoint=s" => \$timepoint,
	    "antsct-base-dir=s" => \$antsCTBaseDir,
	    "antslongct-base-dir=s" => sub {
		my ($opt_name, $opt_value) = @_; 
		$antsCTBaseDir = $opt_value; 
		if ($longitudinalTarget eq "") {
		    $longitudinalTarget = "sst";
		}
	    },
	    "dt-base-dir=s" => \$dtBaseDir,
	    "longitudinal-target=s" => sub { my ($opt_name, $opt_value) = @_; $longitudinalTarget = lc($opt_value); },
	    "output-base-dir=s" => \$outputBaseDir,
	    "cortical-label-system=s" => \$graphCorticalLabelSystem,
	    "cortical-label-def=s" => \$graphCorticalLabelDef
    )
    or die("Error in command line arguments\n");


# Directory for temporary files that is deleted later
my $tmpDir = "";

my $tmpDirBaseName = "${subject}_${timepoint}_dtConnMatPreproc";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = "/tmp/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = $sysTmpDir . "/${tmpDirBaseName}";
}

# Gets removed later, so check we can create this and if not, exit immediately
mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir (maybe it exists from a previous failed run)\n\t";

# Base dir containing DT stuff for this TP
my $tpDTIDir = "${dtBaseDir}/${subject}/${timepoint}";

my $dt = "${tpDTIDir}/dt/${subject}_${timepoint}_DT.nii.gz";

if (! -f "$dt" ) {
    print "\n  No DT for $subject $timepoint\n";
    exit(1);
}

my $dtMask = "${tpDTIDir}/dt/${subject}_${timepoint}_BrainMask.nii.gz";

my $tpOutputDir = "${outputBaseDir}/${subject}/${timepoint}/connMat_${graphCorticalLabelSystem}";     

if (! -d $tpOutputDir) {
    mkpath($tpOutputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $tpOutputDir\n\t";
}

my $t1Brain = "";
my $t1Mask = "";
my $faT1 = "";
my $jlfLabels = "";
my $graphLabels = "";

# For output in SST space only
my $dtToSSTComposedWarp = "";

my $dtToT1DistCorrRoot = "";

if ($longitudinalTarget eq "sst") {
    ($t1Brain, $t1Mask, $faT1, $jlfLabels, $graphLabels, $dtToSSTComposedWarp) = 
	getTargetSpaceImages($antsCTBaseDir, $dtBaseDir, $subject, $timepoint, $graphCorticalLabelSystem, $longitudinalTarget);
    
    if (! -f $dtToSSTComposedWarp ) {
	die("\n  Could not create SST warp for $subject $timepoint");
    }
}
else {
    ($t1Brain, $t1Mask, $faT1, $jlfLabels, $graphLabels, $dtToT1DistCorrRoot) = 
	getTargetSpaceImages($antsCTBaseDir, $dtBaseDir, $subject, $timepoint, $graphCorticalLabelSystem, $longitudinalTarget);
    
    if (! -f "${dtToT1DistCorrRoot}0GenericAffine.mat" ) {
	die("\n  Missing DT -> T1 distortion correction warp for $subject $timepoint");
    }
}

# Check we got all the images we need
if (! -f $t1Brain ) {
    die("\n  No T1 for $subject $timepoint");
}
if (! -f $faT1 ) {
    die("\n  No FA normalized to T1 for $subject $timepoint");
}
if (! -f $t1Mask ) {
    die("\n  No T1 brain mask for $subject $timepoint");
}
if (! -f $jlfLabels ) {
    die("\n  No T1 JLF labels for $subject $timepoint");
}
if (! -f $graphLabels ) {
    die("\n  No graph labels for $subject $timepoint");
}

# Root for output we will create
my $tpOutputRoot = "${tpOutputDir}/${subject}_${timepoint}_";

my $antsVersion = `${antsPath}antsRegistration --version`;

print "ANTs version: ${antsVersion}\n";

print "\nCreating tracking masks from JLF and FA\n";

my ($corticalMask, $exclusionMask, $corticalMaskEdits) = createTrackingMasks($jlfLabels, $faT1, $tpOutputRoot);

print "\nCreating graph nodes, propagating\n  $graphLabels\nto cortical mask\n";

my ($graphNodes, $corticalLabelDiffMask) = createGraphNodes($graphLabels, $graphCorticalLabelDef, $corticalMask, $tpOutputRoot);

my $labeledTractExtent = createTrackingExtentLabeledMask($graphNodes, $corticalMask, $exclusionMask, $tpOutputRoot);

my $dtToReferenceWarpOption = "";

if ($longitudinalTarget eq "sst") {
    $dtToReferenceWarpOption = "--composed-warp $dtToSSTComposedWarp"
}
else {
    $dtToReferenceWarpOption = "--dist-corr-warp-root $dtToT1DistCorrRoot";
}

# Set paths here so they can't get altered by user profiles
my $connmatScriptCmd = qq{
${Bin}/dtConnMat.pl \\
    --dt $dt \\
    --mask $dtMask \\
    --reference-image $t1Brain \\
    ${dtToReferenceWarpOption} \\
    --exclusion-image $exclusionMask \\
    --label-image $graphNodes \\
    --label-def $graphCorticalLabelDef \\
    --output-root ${tpOutputRoot} \\
    --seed-spacing 1 \\
    --seed-fa-thresh 0.2 \\
    --curve-thresh 80 \\
    --compute-scalars 1 \\
    
};

print "\n--- dtConnMat.pl call --- \n${connmatScriptCmd}\n\n";

system("$connmatScriptCmd");

system("rm -f $tmpDir/*");
system("rmdir $tmpDir");


#
# Get images that we need to compute the connectivity in the target space.
#
# my ($t1Brain, $t1Mask, $faT1, $jlfLabels, $graphLabels, $dtToTargetWarp) = getTargetSpaceImages($antsCTBaseDir, $dtBaseDir, $subject, $timepoint, $graphCorticalLabelSystem, $longitudinalTarget)
#
# For longitudinal processing in the SST space, $dtToTargetWarp will be an image, otherwise it will be a string containing
# the distortion correction warp root.
#
# This doesn't check all images for correctness but will return an empty array if the T1 data cannot be found
#
sub getTargetSpaceImages {

    my ($antsCTBaseDir, $dtBaseDir, $subject, $timepoint, $graphCorticalLabelSystem, $longitudinalTarget) = @_;
    
    my $t1Brain = "";
    my $t1Mask = "";
    my $jlfLabels = "";
    my $graphLabels = "";
    my $dtToTargetWarp = "";
    # A session FA warped to T1, space, or for SST output, the average of all time point FA
    my $faT1 = "";

    my $tpDTIDir = "${dtBaseDir}/${subject}/${timepoint}";

    if ($longitudinalTarget) {
	
	my $tpAntsCTDir = `ls -d ${antsCTBaseDir}/${subject}/* | grep "_${timepoint}_"`;
	
	chomp($tpAntsCTDir);

	my $tpSeg = `ls ${tpAntsCTDir} | grep "BrainSegmentation.nii.gz"`;
	
	chomp($tpSeg);
	
	if (! -f "${tpAntsCTDir}/$tpSeg") {
	    print "\n  Incomplete or missing ANTsLongCT data for $subject $timepoint \n";
	    return;
	}
	
	$tpSeg =~ m/(${subject}_${timepoint}_.*)BrainSegmentation.nii.gz/;
	
	my $tpAntsCTFileRoot = $1;
	
	# Full root to antsLongCT output for this time point
	my $tpAntsCTOutputRoot = "${tpAntsCTDir}/${tpAntsCTFileRoot}";
	
	if ($longitudinalTarget eq "sst") {
	    
	    my $sstOutputRoot = "${antsCTBaseDir}/${subject}/${subject}_SingleSubjectTemplate/T_template";
	    
	    $t1Brain = "${sstOutputRoot}ExtractedBrain0N4.nii.gz";

	    $t1Mask = "${sstOutputRoot}BrainExtractionMask.nii.gz";
	    
	    $jlfLabels = "${sstOutputRoot}PG_antsLabelFusionLabels.nii.gz";

	    $graphLabels = "${sstOutputRoot}${graphCorticalLabelSystem}.nii.gz";
	    
	    my $t1ToSSTWarp = "${tpAntsCTOutputRoot}TemplateToSubject0Warp.nii.gz";
	    my $t1ToSSTAffine = "${tpAntsCTOutputRoot}TemplateToSubject1GenericAffine.mat";
	    
	    my $distCorrInvWarp = "${tpDTIDir}/distCorr/${subject}_${timepoint}_DistCorr1InverseWarp.nii.gz";
	    my $distCorrAffine = "${tpDTIDir}/distCorr/${subject}_${timepoint}_DistCorr0GenericAffine.mat";
	    
	    $dtToTargetWarp = "${tmpDir}/${subject}_${timepoint}_tractWarpToSST.nii.gz";

	    my $dtMask = "${tpDTIDir}/dt/${subject}_${timepoint}_BrainMask.nii.gz";

	    system("${antsPath}antsApplyTransforms -d 3 -i $t1Brain -r $dtMask -t [${distCorrAffine}, 1] -t $distCorrInvWarp -t $t1ToSSTAffine -t $t1ToSSTWarp -o [${dtToTargetWarp}, 1] --verbose --float");

	    # subject average FA in SST space
	    $faT1 = "${dtBaseDir}/${subject}/singleSubjectTemplate/${subject}_AverageFA.nii.gz";

	}
	elsif ($longitudinalTarget eq "session") {
	    
	    $t1Brain = "${tpAntsCTOutputRoot}ExtractedBrain0N4.nii.gz";
	    
	    $t1Mask = "${tpAntsCTOutputRoot}BrainExtractionMask.nii.gz";
	    
	    $jlfLabels = "${tpAntsCTOutputRoot}PG_antsLabelFusionLabels.nii.gz";

	    $graphLabels = "${tpAntsCTOutputRoot}${graphCorticalLabelSystem}.nii.gz";
	    
	    $dtToTargetWarp = "${tpDTIDir}/distCorr/${subject}_${timepoint}_DistCorr";

	    $faT1 = "${tpDTIDir}/dtNorm/${subject}_${timepoint}_FANormalizedToStructural.nii.gz";

	}
	else {
	    die("\n  Unrecognized longitudinal target image $longitudinalTarget\n");
	}
    }
    else {
	
	# Cross-sectional pipeline, target is T1
	
	my $tpAntsCTOutputRoot = "${antsCTBaseDir}/${subject}/${timepoint}/${subject}_${timepoint}_";
	
	$t1Brain = "${tpAntsCTOutputRoot}ExtractedBrain0N4.nii.gz";

	# Check for correct input
	if (! -f "$t1Brain") {
	    print "\n  Incomplete or missing ANTsCT data for $subject $timepoint \n";
	    return;
	}
	
	$t1Mask = "${tpAntsCTOutputRoot}BrainExtractionMask.nii.gz";
	
	$jlfLabels = "${tpAntsCTOutputRoot}PG_antsLabelFusionLabels.nii.gz";

	$graphLabels = "${tpAntsCTOutputRoot}${graphCorticalLabelSystem}.nii.gz";

	$dtToTargetWarp = "${tpDTIDir}/distCorr/${subject}_${timepoint}_DistCorr";

	$faT1 = "${tpDTIDir}/dtNorm/${subject}_${timepoint}_FANormalizedToStructural.nii.gz";
	
    }

    return ($t1Brain, $t1Mask, $faT1, $jlfLabels, $graphLabels, $dtToTargetWarp);
    
}


#
# Creates images $corticalMask and $wmMask from the JLF labels and FA in T1w space. 
# The $wmMask contains the WM labels plus connected areas of FA > 0.25 that border cortex. 
#
# The cortical mask is cortical labeled voxels after the extra WM is added. 
#
# The additional WM is designed to prevent false positive connections with when WM was mislabeled as cortex, 
# but we have to be careful to avoid introducing holes into the cortex, which will cause false negatives. 
# Therefore the additional WM is constrained to not reach the edge of the cortical mask.
#
# The first two args are the JLF, labels and FA in T1w space, the rest are file names of 
# images to be created.
#
# The exclusion mask stops tracts that go more than one voxel outside of the WM mask. This allows the 
# streamlines to reach cortical targets, but also lets them go one voxel (in T1w) space outside the WM mask,
# without being terminated.
#
# my ($corticalMask, $exclusionMask, $corticalMaskEdits) = createTrackingMasks($jlfLabels, $faT1, $outputRoot)
#
#
sub createTrackingMasks {

    my ($jlfLabels, $faT1, $outputRoot) = @_;

    # Images to be returned
    my $corticalMask = "${outputRoot}CorticalMask.nii.gz";
    my $exclusionMask = "${outputRoot}ExclusionMask.nii.gz";
    my $corticalMaskEdits = "${outputRoot}CorticalMaskEdits.nii.gz";

    my $tmpOutputRoot = "${tmpDir}/trackingMasks_";

    system("conmat -outputroot ${tmpOutputRoot}conmatCortex_ -targetfile $jlfLabels -targetnamefile $jlfCorticalLabelDef -outputnodes");

    my $corticalLabelMask = "${tmpOutputRoot}corticalJLFMask.nii.gz";
    
    system("${antsPath}ThresholdImage 3 ${tmpOutputRoot}conmatCortex_nodes.nii.gz $corticalLabelMask 1 Inf");
    
    system("conmat -outputroot ${tmpOutputRoot}conmatWM_ -targetfile $jlfLabels -targetnamefile $jlfWMLabelDef -outputnodes");
    
    my $wmLabelMask = "${tmpOutputRoot}wmJLFMask.nii.gz";

    system("${antsPath}ThresholdImage 3 ${tmpOutputRoot}conmatWM_nodes.nii.gz $wmLabelMask 1 Inf");

    my $faMask = "${tmpOutputRoot}faMask.nii.gz";

    system("${antsPath}ThresholdImage 3 $faT1 $faMask 0.25 Inf");

    # Want to avoid adding so much extra to the WM mask that we create holes in the cortical labels
    # One problem is patches of junk FA around the edge of the brain, take large connected components only
    #
    system("${antsPath}LabelClustersUniquely 3 $faMask $faMask 10000");
    system("${antsPath}ThresholdImage 3 $faMask $faMask 1 Inf");
  
    # Further constrain additional WM, add FA > 0.25 only if adjacent to JLF cerebral WM
    my $wmLabelMaskDilated = "${tmpOutputRoot}wmLabelMaskDilated.nii.gz";

    # This dilation sets a global upper limit the maximum extent of WM addition
    system("${antsPath}ImageMath 3 $wmLabelMaskDilated MD $wmLabelMask 2");
    system("${antsPath}ImageMath 3 $faMask m $faMask $wmLabelMaskDilated");

    # Don't allow holes in cortex, so extent of WM addition is constrained to
    # the erosion of GM + WM.
    #
    # This also removes dilation into subcortical structures, eg thalamus or ventricles. 
    #
    # This ought to prevent holes but might fail in some cases, eg it might allow holes in
    # medial GM if the other hemisphere label is close enough and FA is above the threshold in the cortex 
    my $wmAndGMErode = "${tmpOutputRoot}wmAndGMErode.nii.gz";
    system("${antsPath}ImageMath 3 ${tmpOutputRoot}wmAndGM.nii.gz + ${wmLabelMask} ${corticalLabelMask}");
    system("${antsPath}ImageMath 3 $wmAndGMErode ME ${tmpOutputRoot}wmAndGM.nii.gz 1");

    system("${antsPath}ImageMath 3 $faMask m $faMask $wmAndGMErode");

    # Surviving WM (FA > 0.25) that overlaps cortex
    system("${antsPath}ImageMath 3 $faMask m $faMask $corticalLabelMask");

    system("${antsPath}ImageMath 3 ${tmpOutputRoot}faLabelPlusThresh.nii.gz + $faMask $wmLabelMask");

    # Remove anything not connected to cerebral WM 
    system("${antsPath}LabelClustersUniquely 3 ${tmpOutputRoot}faLabelPlusThresh.nii.gz ${tmpOutputRoot}clusters.nii.gz 20000");

    # final WM mask
    my $wmMask = "${tmpOutputRoot}wmMask.nii.gz";
    
    system("${antsPath}ThresholdImage 3 ${tmpOutputRoot}clusters.nii.gz $wmMask 1 Inf");

    if (! -f $wmMask ) {
	die("\n  Could not create WM mask $wmMask");
    }

    my $corticalMaskTmp = "${tmpOutputRoot}corticalMaskTmp.nii.gz";

    # Update the cortical mask
    #
    # We'll do a last check for coverage before making this the final cortical mask
    system("${antsPath}ImageMath 3 $corticalMaskTmp - $corticalLabelMask $wmMask");
    system("${antsPath}ThresholdImage 3 $corticalMaskTmp $corticalMaskTmp 1 1");

    # Now make a WM tract termination mask
    
    my $wmMaskDilated = "${tmpOutputRoot}wmMaskDilated.nii.gz";
    
    system("${antsPath}ImageMath 3 $wmMaskDilated MD $wmMask 1");
    
    # Termination mask for Camino is the inversion of the dilated mask
    system("${antsPath}ThresholdImage 3 $wmMaskDilated $exclusionMask 1 1 0 1");
    
    if (! -f $exclusionMask ) {
	die("\n  Could not create exclusion mask $exclusionMask");
    }

    # Now check that the boundary of the termination mask has cortical labels in voxels that are next to cortex
    my $trackingEdgeMask = "${tmpOutputRoot}trackingBoundaryMask.nii.gz";
    
    system("${antsPath}ImageMath 3 $trackingEdgeMask - $wmMaskDilated $wmMask");

    my $corticalMaskDilated = "${tmpOutputRoot}corticalMaskDilated.nii.gz";

    system("${antsPath}ImageMath 3 $corticalMaskDilated MD $corticalMaskTmp 1");

    # The intersection of these two should contain the cortical labels
    my $possibleHoles = "${tmpOutputRoot}fillHolesMask.nii.gz";

    system("${antsPath}ImageMath 3 $possibleHoles m $corticalMaskDilated $trackingEdgeMask");

    # Add these to the cortical mask
    system("${antsPath}ImageMath 3 $corticalMask addtozero $corticalMaskTmp $possibleHoles");

    if (! -f $corticalMask ) {
	die("\n  Could not create cortical mask $corticalMask");
    }

    # For QC purposes, record edits made to the cortical mask
    # 1 = in JLF mask, but not in final mask
    # 2 = in final mask, but not in JLF mask
    # 3 = in both
    system("${antsPath}ImageMath 3 $corticalMaskEdits m $corticalMask 2");
    system("${antsPath}ImageMath 3 $corticalMaskEdits + $corticalMaskEdits $corticalLabelMask");

    if (! -f $corticalMaskEdits ) {
	die("\n  Could not create QC image $corticalMaskEdits");
    }

    return ($corticalMask, $exclusionMask, $corticalMaskEdits);

}


#
# Create $graphNodes for some label set. 
#
# Also make a QC image with the overlap between the graph node mask and the input cortical mask.
#
# my ($graphNodes, $corticalLabelDiffMask) = createGraphNodes($labelImage, $corticalLabelDef, $corticalMask, $outputRoot)
# 
sub createGraphNodes {
    
    my ($labelImage, $corticalLabelDef, $corticalMask, $outputRoot) = @_;

    my $graphNodes = "${outputRoot}GraphNodes.nii.gz";

    my $tmpOutputRoot = "${tmpDir}/createGraphNodes_";

    system("conmat -outputroot ${tmpOutputRoot}conmat_ -targetfile $labelImage -targetnamefile $corticalLabelDef -outputnodes");
    
    # Use this to compare to final version
    my $corticalLabelsOrig = "${tmpOutputRoot}conmat_nodes.nii.gz";

    my $corticalMaskOrig = "${tmpOutputRoot}corticalMaskOrig.nii.gz";

    system("${antsPath}ThresholdImage 3 $corticalLabelsOrig $corticalMaskOrig 1 Inf");

    my $corticalLabelsTmp = "${tmpOutputRoot}corticalLabelsTmp.nii.gz";
    
    # Propagate the labels to the mask (fills unlabeled voxels in mask)
    #
    # Last two parameters add stopping criteria and topology constraints
    system("${antsPath}ImageMath 3 $corticalLabelsTmp PropagateLabelsThroughMask $corticalMask $corticalLabelsOrig 100 1");

    # Don't need to separately mask the labels to remove labeled voxels outside the mask, PropagateLabels does this

    system("cp $corticalLabelsTmp $graphNodes");

    if (! -f $graphNodes ) {
	die("\n  Could not create graph nodes $graphNodes");
    }
    
    my $corticalLabelDiffMask = "${outputRoot}CorticalLabelSystemDiffMask.nii.gz";

    # QC - show overlap of input label image and graph nodes. Large differences may indicate poor registration
    #
    # 1 = in original mask, but not graph nodes
    # 2 = in graph nodes, but not original mask
    # 3 = in both
    system("${antsPath}ImageMath 3 $corticalLabelDiffMask m $corticalMask 2");
    system("${antsPath}ImageMath 3 $corticalLabelDiffMask + $corticalLabelDiffMask $corticalMaskOrig");

    return ($graphNodes, $corticalLabelDiffMask);
    
}


#
# For QC
#
# Create $labeledTractExtent, this image is > 1 in all voxels that tracts can reach.
#
# This can be visualized to see that tracts can reach the graph nodes
#
# my $labeledTractExtent = createTrackingExtentLabeledMask($graphNodes, $corticalMask, $exclusionMask, $outputRoot)
#
sub createTrackingExtentLabeledMask {

    my ($graphNodes, $corticalMask, $exclusionMask, $outputRoot) = @_;

    my $labeledTractExtent = "${outputRoot}LabeledTractInclusionMask.nii.gz";
    
    my $tmpOutputRoot = "${tmpDir}/createTrackingExtentLabeledMask_";

    my $tractExtentMask = "${tmpOutputRoot}TractExtentMask.nii.gz";

    system("${antsPath}ThresholdImage 3 $exclusionMask $tractExtentMask 0 0");

    my $corticalBoundaryVoxels = "${tmpOutputRoot}corticalBoundaryVoxels.nii.gz";

    system("${antsPath}ImageMath 3 $corticalBoundaryVoxels m $tractExtentMask $graphNodes");

    # Find the maximum cortical label; WM label is max + 1

    my $minMaxMeanOutput = `MeasureMinMaxMean 3 $corticalBoundaryVoxels`; 
    
    $minMaxMeanOutput =~ m/Max : \[(\d+)\]/;

    my $maxLabel = $1 + 1;

    my $wmLabeled = "${tmpOutputRoot}wmLabeled.nii.gz";

    system("${antsPath}ImageMath 3 $wmLabeled m $tractExtentMask $maxLabel");

    system("${antsPath}ImageMath 3 $labeledTractExtent addtozero $corticalBoundaryVoxels $wmLabeled ");

    if (! -f $labeledTractExtent ) {
	die("\n  Could not create QC image $labeledTractExtent");
    }

    return $labeledTractExtent;
}
