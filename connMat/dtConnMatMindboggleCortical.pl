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

# CSV file containing Mindboggle label definitions
# These may be cortical or cortical + subcortial, but are assumed to be mindboggle
# This script is not general to different label schemes because it uses label > 100 
# as a cortical mask
my $corticalLabelDef = "${Bin}/mindBoggleCorticalGraphNodes.csv";

# Used to make an exclusion ROI
my $wmLabelDef = "${Bin}/mindBoggleWMLabels.csv";


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
     Base output directory. Output is organized under this as subj/tp/connMat.


  Options:

   --longitudinal-target
     Required only for longitudinal data. Either "session" to evaluate connectivity in the intra-session 
     T1, or "sst", to do it in the SST space. If "session", graph nodes must be defined for each time point.
  


  Some pre-processing is done at run time so it is best to run this script from a qlogin session, or qsub it. 
  You may need to set CAMINO_HEAP_SIZE first to allocate enough RAM for the creation of the nodes.

  The target T1 image / SST should be labeled using the pseudo-geodesic JLF.

  DT data is read from

    \${dt-base-dir}/subj/tp/

  The DTI preprocessing pipeline should have been run, such that dt/ and distCorr/ exist; 
  these will be used to do the tracking and transfer the results to T1 space.

  
  Output:

    Prefixed with outputBaseDir/subj/tp/connMat/subj_tp_ :

      GraphNodes.nii.gz     
 
    Mindboggle cortical nodes, masked to exclude high FA without introducing holes.


      ExclusionMask.nii.gz

    A mask defined by the dilation of the JLF WM segmentation and adjoining regions of high FA.
    This constrains tracking to white matter, but still allows lines to reach the nodes.

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
	    "output-base-dir=s" => \$outputBaseDir
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

my $tpOutputDir = "${outputBaseDir}/${subject}/${timepoint}/connMat";     

if (! -d $tpOutputDir) {
    mkpath($tpOutputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $tpOutputDir\n\t";
}

my $t1Brain = "";
my $t1Mask = "";
my $faT1 = "";
my $jlfLabels = "";

# For output in SST space only
my $dtToSSTComposedWarp = "";

my $dtToT1DistCorrRoot = "";

if ($longitudinalTarget eq "sst") {
    ($t1Brain, $t1Mask, $faT1, $jlfLabels, $dtToSSTComposedWarp) = 
	getTargetSpaceImages($antsCTBaseDir, $dtBaseDir, $subject, $timepoint, $longitudinalTarget);
    
    if (! -f $dtToSSTComposedWarp ) {
	die("\n  Could not create SST warp for $subject $timepoint");
    }
}
else {
    ($t1Brain, $t1Mask, $faT1, $jlfLabels, $dtToT1DistCorrRoot) = 
	getTargetSpaceImages($antsCTBaseDir, $dtBaseDir, $subject, $timepoint, $longitudinalTarget);
    
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

# Root for output we will create
my $tpOutputRoot = "${tpOutputDir}/${subject}_${timepoint}_";

my $wmMask = "${tmpDir}/${subject}_${timepoint}_wmMask.nii.gz";

createWMMask($jlfLabels, $corticalLabelDef, $wmLabelDef, $faT1, $wmMask);

# Exclusion ROI is the negation of the dilated WM mask
my $exclusionMask = "${tpOutputDir}/${subject}_${timepoint}_ExclusionMask.nii.gz";

createExclusionMask($wmMask, $exclusionMask);

my $graphNodes = "${tpOutputDir}/${subject}_${timepoint}_GraphNodes.nii.gz";

createGraphNodes($jlfLabels, $corticalLabelDef, $wmMask, $graphNodes);

my $scriptToRun = "${tpOutputDir}/connMat_${subject}_${timepoint}.sh";

my $antsVersion = `${antsPath}antsRegistration --version`;

my $dtToReferenceWarpOption = "";

if ($longitudinalTarget eq "sst") {
    $dtToReferenceWarpOption = "--composed-warp $dtToSSTComposedWarp"
}
else {
    $dtToReferenceWarpOption = "--dist-corr-warp-root $dtToT1DistCorrRoot";
}

print "ANTs version: ${antsVersion}\n";

# Set paths here so they can't get altered by user profiles
my $connmatScriptCmd = qq{
${Bin}/dtConnMat.pl \\
    --dt $dt \\
    --mask $dtMask \\
    --reference-image $t1Brain \\
    ${dtToReferenceWarpOption} \\
    --exclusion-image $exclusionMask \\
    --label-image $graphNodes \\
    --label-def $corticalLabelDef \\
    --output-root ${tpOutputRoot} \\
    --seed-spacing 1 \\
    --seed-fa-thresh 0.2 \\
    --curve-thresh 80 \\
    --compute-scalars 1 \\
    
};

print "--- dtConnMat.pl call --- \n${connmatScriptCmd}\n\n";

system("$connmatScriptCmd");

system("rm -f $tmpDir/*");
system("rmdir $tmpDir");


#
# Get images that we need to compute the connectivity in the target space.
#
# my ($t1Brain, $t1Mask, $faT1, $jlfLabels, $dtToTargetWarp) = getTargetSpaceImages($antsCTBaseDir, $dtBaseDir, $subject, $timepoint, $longitudinalTarget)
#
# For longitudinal processing in the SST space, $dtToTargetWarp will be an image, otherwise it will be a string containing
# the distortion correction warp root.
#
# This doesn't check all images for correctness but will return an empty array if the T1 data cannot be found
#
sub getTargetSpaceImages {

    my ($antsCTBaseDir, $dtBaseDir, $subject, $timepoint, $longitudinalTarget) = @_;
    
    my $t1Brain = "";
    my $t1Mask = "";
    my $jlfLabels = "";
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

	$dtToTargetWarp = "${tpDTIDir}/distCorr/${subject}_${timepoint}_DistCorr";

	$faT1 = "${tpDTIDir}/dtNorm/${subject}_${timepoint}_FANormalizedToStructural.nii.gz";
	
    }

    return ($t1Brain, $t1Mask, $faT1, $jlfLabels, $dtToTargetWarp);
    
}


#
# creates image $wmMask containing the WM labels plus connected areas of FA > 0.25 
# that border cortex. 
#
# This is used to constrain tracking and to update the cortical nodes by excluding 
# high FA voxels.
#
# createExclusionMask($jlfLabels, $corticalLabelDef, $wmLabelDef, $wmMask)
#
sub createWMMask  {

    my ($jlfLabels, $corticalLabelDef, $wmLabelDef, $faT1, $wmMask) = @_;

    my $tmpOutputRoot = "${tmpDir}/" . fileparse($wmMask, (".nii", ".nii.gz")) . "_";

    system("conmat -outputroot ${tmpOutputRoot}conmatCortex_ -targetfile $jlfLabels -targetnamefile $corticalLabelDef -outputnodes");

    my $corticalNodeMask = "${tmpOutputRoot}corticalNodeMask.nii.gz";
    
    system("${antsPath}ThresholdImage 3 ${tmpOutputRoot}conmatCortex_nodes.nii.gz $corticalNodeMask 1 Inf");
    
    system("conmat -outputroot ${tmpOutputRoot}conmatWM_ -targetfile $jlfLabels -targetnamefile $wmLabelDef -outputnodes");
    
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
    system("${antsPath}ImageMath 3 $wmLabelMaskDilated MD $wmLabelMask 2");
    system("${antsPath}ImageMath 3 $faMask m $faMask $wmLabelMaskDilated");

    # Don't allow holes in cortex, maximum extent of WM addition is constrained to
    # the erosion of GM + WM.
    #
    # This by itself ought to be good but might fail in some cases, eg it might allow holes in
    # medial GM if the other hemisphere label is close enough
    my $wmAndGMErode = "${tmpOutputRoot}wmAndGMErode.nii.gz";
    system("${antsPath}ImageMath 3 ${tmpOutputRoot}wmAndGM.nii.gz + ${wmLabelMask} ${corticalNodeMask}");
    system("${antsPath}ImageMath 3 $wmAndGMErode ME ${tmpOutputRoot}wmAndGM.nii.gz 1");

    system("${antsPath}ImageMath 3 $faMask m $faMask $wmAndGMErode");

    # Surviving WM (FA > 0.25) that overlaps cortex. Remove cortical labels from such voxels
    system("${antsPath}ImageMath 3 $faMask m $faMask $corticalNodeMask");

    system("${antsPath}ImageMath 3 ${tmpOutputRoot}faLabelPlusThresh.nii.gz + $faMask $wmLabelMask");

    system("${antsPath}LabelClustersUniquely 3 ${tmpOutputRoot}faLabelPlusThresh.nii.gz ${tmpOutputRoot}clusters.nii.gz 20000");

    system("${antsPath}ThresholdImage 3 ${tmpOutputRoot}clusters.nii.gz $wmMask 1 Inf");

    if (! -f $wmMask ) {
	die("\n  Could not create WM mask $wmMask");
    }

}


#
# Cortical label voxels inside the WM mask will be removed. 
#
# createGraphNodes($jlfLabels, $corticalLabelDef, $wmMask, $graphNodes)
#
sub createGraphNodes {

    my ($jlfLabels, $corticalLabelDef, $wmMask, $graphNodes) = @_;

    my $tmpOutputRoot = "${tmpDir}/" . fileparse($graphNodes, (".nii", ".nii.gz")) . "_";
    
    system("conmat -outputroot ${tmpOutputRoot}conmat_ -targetfile $jlfLabels -targetnamefile $corticalLabelDef -outputnodes");

    # This mask is created in several steps, the end result being that it's 1 if the voxel has a cortical label but is probably
    # WM. It is used to remove such voxels from the final nodes
    my $corticalNodeMask = "${tmpOutputRoot}corticalNodeMask.nii.gz";
    my $corticalNodeExclusionMask = "${tmpOutputRoot}corticalExclusionNodeMask.nii.gz";
 
    system("${antsPath}ThresholdImage 3 ${tmpOutputRoot}conmat_nodes.nii.gz $corticalNodeMask 1 Inf");
 
    system("${antsPath}ImageMath 3 $corticalNodeMask m $wmMask $corticalNodeMask");

    system("${antsPath}ThresholdImage 3 $corticalNodeMask $corticalNodeExclusionMask 1 1 0 1");

    system("${antsPath}ImageMath 3 $graphNodes m ${tmpOutputRoot}conmat_nodes.nii.gz $corticalNodeExclusionMask");

    if (! -f $graphNodes ) {
	die("\n  Could not create nodes $graphNodes");
    }

}



#
# Makes an exclusion (tract termination) mask from the negation of (dilated JLF WM + the nodes we are interested in).
#
# Streamlines will be able to propagate 
# 
# createExclusionMask($wmMask, $exclusionMask)
#
sub createExclusionMask  {

    my ($wmMask, $exclusionMask) = @_;

    my $tmpOutputRoot = "${tmpDir}/" . fileparse($exclusionMask, (".nii", ".nii.gz")) . "_";
    
    my $tmpMask = "${tmpOutputRoot}wmMaskDilated.nii.gz";
    
    system("${antsPath}ImageMath 3 $tmpMask MD $wmMask 1");
    
    system("${antsPath}ThresholdImage 3 $tmpMask $exclusionMask 1 1 0 1");
    
    if (! -f $exclusionMask ) {
	die("\n  Could not create tracking exclusion mask $exclusionMask");
    }
}
