#!/usr/bin/perl -w
#
# Wrapper script for processing DTI for a particular scan. Runs all the processing steps in series.
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
      --antsct-output-root 
     
      [ options ]

  This script is for processing non-topup data, using a T1 that has been processed with antsCorticalThickness, a local 
  template, and optionally a standard template.

  Required args:

   --input-dir
     Input directory containing the DWI data, bvals, and bvecs.

   --output-dir
     Output directory. Subdirectories with the various processing steps will be created (see below).

   --output-file-root
     Prepended onto output files. A sensible value of this would be \${subject}_\${tp}.

   --antsct-output-root
     Path to output from antsCT, such that \${root}ExtractedBrain0N4.nii.gz exists. Other antsCT output is also required,
     including the brain mask and the warps to the local template space.


  Options:

   --acq-params
     Acquition parameters file, required for eddy. The index file will be auto-generated, so it is assumed that the 
     DWI data has a consistent phase encode direction and echo spacing, but the number of measurements 
     is inferred from the bvals.

   --eddy-correct-method
     Either "ANTs" or "FSL". Default ANTs, which does affine registration to a b=0 reference image. With FSL, 
     uses eddy.

   --template
     Template image used as a reference space. The warps from the structural image to this template should exist in the 
     antsCT directory.

   --template-mask
     Mask image in the template space, used to mask deformed DT.

   --standard-template
     Standard template space, used as a reference space. Warps go via the local template.

   --standard-template-mask
     Mask image in the standard space, used to mask deformed DT.

   --standard-template-warp-root 
     Warp root such that root0GenericAffine.mat exists, and optionally root1Warp.nii.gz exists also. The warp should map 
     from the local template to the standard template space.
     

  Output:

   Output is organized in subdirectories under the main output dir (in order of processing):

     merged/ - Merged DWI data. If there are repeats, they get combined into a single file, with 
               associated bvals and bvecs.

     distCorr/ - Distortion correction mapping to T1, and the brain mask for the DTI 
                 imported from the T1.
      
     eddy/ - eddy-corrected data and bvectors.

     dtfit/ - tensors and associated information.

     dtNormalized/ - tensors and associated information normalized to other image spaces: T1, local template, and MNI.

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
my ($antsCTRoot, $inputDir, $outputFileRoot, $outputScanDir);

# Options
my $acqParams = ""; 
my $eddyCorrect = "ants";
my $template = "";
my $templateMask = "";
my $standardTemplate = "";
my $standardTemplateWarpRoot = "";
my $standardTemplateMask = "";

GetOptions ("input-dir=s" => \$inputDir,
	    "output-dir=s" => \$outputScanDir,
	    "output-file-root=s" => \$outputFileRoot,
	    "acq-params=s" => \$acqParams,
	    "eddy-correct-method" => $eddyCorrect,
	    "antsct-output-root=s" => \$antsCTRoot,
	    "template=s" => \$template,
	    "template-mask=s" => \$templateMask,
	    "standard-template=s" => \$standardTemplate,
	    "standard-template-warp-root=s" => \$standardTemplateWarpRoot,
	    "standard-template-mask=s" => \$standardTemplateMask,
    )
    or die("Error in command line arguments\n");

$eddyCorrect = lc($eddyCorrect);

# We will not quit here if output exists, rather check each output subdir and skip those that exist. This allows 
# steps to be redone if needed
if (! -d $outputScanDir ) { 
  mkpath($outputScanDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputScanDir\n\t";
}

# Merge DWI data into single series

print "\nMerge DWI data\n";

my $mergeDataDir = "${outputScanDir}/merged";

my $mergeDataOutputRoot = "${mergeDataDir}/${outputFileRoot}DWI";

if ( -d $mergeDataDir ) {
    print "\n  Merged output directory $mergeDataDir exists, will not merge raw data\n"; 
}
else {
    my $cmd = "${Bin}/mergeDiffusionData.pl --input-dir ${inputDir} --output-root ${mergeDataOutputRoot}";
    
    print "\n--- Merge data command ---\n${cmd}\n---\n";
    
    system("${cmd}") == 0 or die("\nDWI data merge failed\n");
    
}


# Distortion correction to T1

my $distCorrDir = "${outputScanDir}/distCorr";

my $distCorrOutputRoot = "${distCorrDir}/${outputFileRoot}DistCorr";

# produced by distcorr script, used in other functions
my $dwiBrainMask = "${distCorrOutputRoot}DWIBrainMask.nii.gz";

my $t1 = "${antsCTRoot}ExtractedBrain0N4.nii.gz";
my $t1BrainMask = "${antsCTRoot}BrainExtractionMask.nii.gz";

my $mocoB0 = "${distCorrDir}/${outputFileRoot}MocoB0.nii.gz";

if ( -d $distCorrDir ) {
    print "\n  DistCorr output directory $distCorrDir exists, will not correct data\n"; 
}
else {
    my $cmd = "${Bin}/mocoB0.pl --dwi ${mergeDataOutputRoot}.nii.gz --bvecs ${mergeDataOutputRoot}.bvec --bvals ${mergeDataOutputRoot}.bval --output-file $mocoB0";

    print "\n--- Compute B0 for distortion correction ---\n${cmd}\n---\n";

    system("${cmd}");

    $cmd = "${Bin}/distCorrDiffToStructural.pl --structural $t1 --mask $t1BrainMask --dwi $mocoB0 --output-root $distCorrOutputRoot";

    print "\n--- Distortion correction to T1 ---\n${cmd}\n---\n";

    system("${cmd}") == 0 or die("\nDT to structural distortion correction failed\n");

}


# eddy and motion correct

my $eddyDir = "${outputScanDir}/eddy";

my $eddyOutputRoot = "${eddyDir}/${outputFileRoot}DWICorrected";

if ( -d $eddyDir ) {
    print "\n  eddy output directory $eddyDir exists, will not run eddy\n"; 
}
else {

    if ($eddyCorrect eq "ants") {
	my $cmd = "${Bin}/antsEddyCorrect.pl --dwi ${mergeDataOutputRoot}.nii.gz --bvals ${mergeDataOutputRoot}.bval --bvecs ${mergeDataOutputRoot}.bvec --brain-mask $dwiBrainMask --ref $mocoB0 --output-root $eddyOutputRoot";

	print "\n--- eddy ---\n${cmd}\n---\n";
	
	system("${cmd}") == 0 or die("\neddy failed\n");

    }
    elsif ($eddyCorrect eq "fsl") {
	my $cmd = "${Bin}/runEddy.sh -i ${mergeDataOutputRoot}.nii.gz -b ${mergeDataOutputRoot}.bval -r ${mergeDataOutputRoot}.bvec -m $dwiBrainMask -a $acqParams -o $eddyOutputRoot -q 0 -t 0";

	print "\n--- eddy ---\n${cmd}\n---\n";
	
	system("${cmd}") == 0 or die("\neddy failed\n");

	system("mv ${eddyOutputRoot}.eddy_rotated_bvecs  ${eddyOutputRoot}.bvec");
    }
    else {
	die("Unrecognized eddy correction option $eddyCorrect\n\t");
    }
}


# fit the DT

my $dtDir = "${outputScanDir}/dt";

my $dtOutputRoot = "${dtDir}/${outputFileRoot}";

if ( -d $dtDir ) {
    print "\n  DT output directory $dtDir exists, will not run DT fit\n"; 
}
else {

    my $cmd = "${Bin}/fitDT.pl --dwi ${eddyOutputRoot}.nii.gz --mask $dwiBrainMask --bvals ${eddyOutputRoot}.bval --bvecs ${eddyOutputRoot}.bvec --output-root $dtOutputRoot --algorithm wdt";

    print "\n--- DT fit ---\n${cmd}\n---\n";

    system("${cmd}") == 0 or die("\nDT fitting failed\n");

}


# Warp DT to T1, and optionally templates

my $templateDTDir = "${outputScanDir}/dtNorm";

my $templateDT_OutputRoot = "${templateDTDir}/${outputFileRoot}";

if ( -d $templateDTDir ) {
    print "\n output directory $templateDTDir exists, will not run DT warp\n"; 
}
else {

    # Get all possible output here. To skip certain outputs, supply the warps but not the reference image. For example, if you don't specify
    # --structural, you won't get FA etc resampled in T1 space.

    my $dwiBrainMaskCorrected = "${dtOutputRoot}BrainMask.nii.gz";
   
    my $cmd = "${Bin}/dtToTemplate.pl --dt ${dtOutputRoot}DT.nii.gz --mask $dwiBrainMaskCorrected --dist-corr-warp-root $distCorrOutputRoot --output-root $templateDT_OutputRoot --structural $t1 --structural-mask $t1BrainMask";

    if (-f $template) {
	$cmd = $cmd . " --template $template --template-warp-root ${antsCTRoot}SubjectToTemplate";

	if (-f $templateMask ) {
	    $cmd = $cmd . " --template-mask $templateMask";
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
