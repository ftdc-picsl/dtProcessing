#!/usr/bin/perl -w
#
# Wrapper script for processing DTI for a particular scan
#

use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;
use File::Basename;
use Getopt::Long;

my $usage = qq{

  $0  
      --input-dir 
      --output-dir 
      --output-file-root 
      --antslongct-timepoint-output-root
      --antslongct-sst-output-root
     
      [ options ]

  This script is for processing non-topup data, using a T1 that has been processed with antsLongitudinalCorticalThickness.

  The data must have been pre-processed in the cross-sectional pipeline first. The longitudinal processing steps are:

  1. Distortion correction to the timepoint T1. This is done already cross-sectionally, but is recomputed in case the
     T1 image has changed for longitudinal processing.

  2. Mask the corrected DWI image, and fit the DT.

  3. Apply the warps to move the DT to the SST and the standard template space.

  Required args:

   --input-dir
     Input directory containing the pre-processed DWI data from the cross-sectional pipeline.

   --output-dir
     Output directory. Subdirectories with the various processing steps will be created (see below).

   --output-file-root
     Prepended onto output files. A sensible value of this would be \${subject}_\${tp}.

   --antslongct-timepoint-output-root
     Path to output from antsLongCT, such that \${root}ExtractedBrain0N4.nii.gz is the T1 image from the same 
     scanning session as the DTI. 

   --antslongct-sst-output-root
     Path to output from antsLongCT, such that \${root}ExtractedBrain0N4.nii.gz is the SST of all time points.


  Options:

   --group-template
     Group space reference image. Required to resample the DT into the local group template space.

   --group-template-mask
     A brain mask in the group template space. Used to set the background tensors to zero after warping.

   --standard-template
     Standard space reference image. Required to resample the DT into the standard template space. Also the 
     the local group template.

   --standard-template-warp-root 
     Forward warps defining a transform from the local group template to the standard space. 

   --standard-template-mask
     A brain mask in the standard template space. Used to set the background tensors to zero after warping.


  Output:

   Output is organized in subdirectories under the main output dir (in order of processing):

     distCorr/ - Distortion correction mapping to T1, and the brain mask for the DTI 
                 imported from the T1.
      
     dt/ - tensors and associated information.

     dtNormalized/ - tensors and associated information normalized to other image spaces: 
                     T1, SST, group template, standard template

   If any of the above output directories exist, that step is skipped, and the script continues. 


  Requires ANTs and Camino

};


if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

# Get the directories containing programs we need
my ($antsPath, $sysTmpDir) = @ENV{'ANTSPATH', 'TMPDIR'};

if (!$antsPath || ! -f "${antsPath}antsRegistration") {
    die("Script requires ANTSPATH\n\t");
}

my $dtfit = `which dtfit`;

chomp $dtfit;

if (! -f $dtfit) {
    die("Script requires Camino\n\t");
}


# Required args
my ($antsLongCT_TP_OutputRoot, $antsLongCT_SST_OutputRoot, $inputDir, $outputFileRoot, $outputScanDir);

# Options

# Group template passed to antsLongCT
# Warp root defined by antsLongCT
my $groupTemplate = "";
my $groupTemplateMask = "";

# MNI template
my $standardTemplate = "";
# Warp from group template to MNI
my $standardTemplateWarpRoot="";
my $standardTemplateMask = "";

GetOptions ("input-dir=s" => \$inputDir,
	    "output-dir=s" => \$outputScanDir,
	    "output-file-root=s" => \$outputFileRoot,
	    "antslongct-timepoint-output-root=s" => \$antsLongCT_TP_OutputRoot,
	    "antslongct-sst-output-root=s" => \$antsLongCT_SST_OutputRoot,
	    "group-template=s" => \$groupTemplate,
	    "group-template-mask=s" => \$groupTemplateMask,
	    "standard-template=s" => \$standardTemplate,
	    "standard-template-warp-root=s" => \$standardTemplateWarpRoot,
	    "standard-template-mask=s" => \$standardTemplateMask,
    )
    or die("Error in command line arguments\n");


my $sst = "${antsLongCT_SST_OutputRoot}ExtractedBrain0N4.nii.gz";
# Warp from T1 to SST
my $sstWarpRoot = "${antsLongCT_TP_OutputRoot}SubjectToTemplate";
my $sstMask = "${antsLongCT_SST_OutputRoot}BrainExtractionMask.nii.gz";

# Warp from SST to group template
my $groupTemplateWarpRoot = "${antsLongCT_SST_OutputRoot}SubjectToTemplate";

# We will not quit here if output exists, rather check each output subdir and skip those that exist. This allows 
# steps to be redone if needed
if (! -d $outputScanDir ) { 
  mkpath($outputScanDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputScanDir\n\t";
}

# Distortion correction to timepoint T1

my $distCorrDir = "${outputScanDir}/distCorr";

my $distCorrOutputRoot = "${distCorrDir}/${outputFileRoot}DistCorr";

# produced by distcorr script, used in other functions
my $dwiBrainMask = "${distCorrOutputRoot}DWIBrainMask.nii.gz";

my $t1 = "${antsLongCT_TP_OutputRoot}ExtractedBrain0N4.nii.gz";
my $t1BrainMask = "${antsLongCT_TP_OutputRoot}BrainExtractionMask.nii.gz";

if ( -d $distCorrDir ) {
    print "\n  DistCorr output directory $distCorrDir exists, will not correct data\n"; 
}
else {

    # From the cross-sectional output
    my $mocoB0 = `ls ${inputDir}/distCorr/*MocoB0.nii.gz`;

    chomp($mocoB0);

    my $cmd = "${Bin}/distCorrDiffToStructural.pl --structural $t1 --mask $t1BrainMask --dwi $mocoB0 --output-root $distCorrOutputRoot";

    print "\n--- Distortion correction to T1 ---\n${cmd}\n---\n";

    system("${cmd}") == 0 or die("\nDT to structural distortion correction failed\n");

}


# fit the DT

my $dtDir = "${outputScanDir}/dt";

my $dtOutputRoot = "${dtDir}/${outputFileRoot}";

if ( -d $dtDir ) {
    print "\n  DT output directory $dtDir exists, will not run DT fit\n"; 
}
else {

    my $dwiCorrected = `ls ${inputDir}/eddy/*DWICorrected.nii.gz`; 
    my $bvals = `ls ${inputDir}/eddy/*.bval`; 
    my $bvecs = `ls ${inputDir}/eddy/*.bvec`; 

    chomp($dwiCorrected, $bvals, $bvecs);

    my $cmd = "${Bin}/fitDT.pl --dwi $dwiCorrected --mask $dwiBrainMask --bvals $bvals --bvecs $bvecs --output-root $dtOutputRoot --algorithm wdt";

    print "\n--- DT fit ---\n${cmd}\n---\n";

    system("${cmd}") == 0 or die("\nDT fitting failed\n");

}


# Warp DT to T1, and optionally templates

my $templateDT_Dir = "${outputScanDir}/dtNorm";

my $templateDT_OutputRoot = "${templateDT_Dir}/${outputFileRoot}";

if ( -d $templateDT_Dir ) {
    print "\n output directory $templateDT_Dir exists, will not run DT warp\n"; 
}
else {

    # Get all possible output here. To skip certain outputs, supply the warps but not the reference image. 
    # For example, if you don't specify --structural, you won't get the DT resampled in T1 space.

    my $dwiBrainMaskCorrected = "${dtOutputRoot}BrainMask.nii.gz";
   
    my $cmd = "${Bin}/dtToTemplateLong.pl --dt ${dtOutputRoot}DT.nii.gz --mask $dwiBrainMaskCorrected --dist-corr-warp-root $distCorrOutputRoot --output-root $templateDT_OutputRoot --structural $t1 --structural-mask $t1BrainMask";

    $cmd = $cmd . " --single-subject-template $sst --single-subject-template-warp-root $sstWarpRoot --single-subject-template-mask $sstMask";

    $cmd = $cmd . " --group-template-warp-root $groupTemplateWarpRoot";

    if (-f $groupTemplate) {
	$cmd = $cmd . " --group-template $groupTemplate";

	if (-f $groupTemplateMask ) {
	    $cmd = $cmd . " --group-template-mask $groupTemplateMask";
	}
    }

    if (-f $standardTemplate) {
	$cmd = $cmd . " --standard-template $standardTemplate --standard-template-warp-root $standardTemplateWarpRoot";

	if (-f $standardTemplateMask ) {
	    $cmd = $cmd . " --standard-template-mask $standardTemplateMask";
	}

    }

    print "\n--- Warp DT to template space(s) ---\n${cmd}\n---\n";

    system("${cmd}") == 0 or die("\nWarp DT to standard space failed\n");

}
