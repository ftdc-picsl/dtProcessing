#!/usr/bin/perl -w
#
# Wrapper script for calling processScanDTI.pl for some or all of a subject's data
#

use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;
use File::Basename;
use Getopt::Long;


my $templateDir = "/data/grossman/pipedream2018/templates/OASIS";

my $template = "${templateDir}/T_template0_BrainCerebellum.nii.gz";

my $templateToMNI_WarpRoot = "${templateDir}/MNI152/T_template0_ToMNI152";

my $mniTemplate ="${templateDir}/MNI152/MNI152_T1_1mm_brain.nii.gz";

# Masks are used to mask deformed DT 
my $templateMask = "${templateDir}/T_template0_BrainCerebellumMask.nii.gz";

my $mniTemplateMask = "${templateDir}/MNI152/MNI152_1mm_brainMask.nii.gz";

my $antsPath = "/data/grossman/pipedream2018/bin/ants/bin/";

my $caminoDir = "/data/grossman/pipedream2018/bin/camino/bin";

my $acqParamsDir = "/data/grossman/pipedream2018/metadata";


# where to find raw NII DWI data, bvecs, bvals
my $dwiBaseDir = "";

# Base dir for output
my $outputBaseDir = "";

# Base dir for T1 antsCT
my $antsCTBaseDir = "";


my $cores=1;
my $ram="4";
my $submitToQueue=1;
my $eddyCorrectMethod = "ants";

my $usage = qq{

  $0  
      --subject
      --dwi-base-dir
      --anstct-base-dir
      --output-base-dir
     
      [ options ]

  Required args:

   --subject
     Subject ID.

   --dwi-base-dir
     Input base directory where subject/timepoint/DWI/ has the raw DWI data and bvecs / bvals.

   --anstct-base-dir
     Base directory for antsCT, where subject/timepoint/ contains cross-sectional antsCT output.

   --output-base-dir
     Output base directory, output is organized by subject / timepoint under this.



  Options:

   --timepoints  
     Timepoint(s) to process. If no time points are provided, it means process all available timepoints.

   --qsub
     Submit processing jobs to qsub (default = $submitToQueue).

   --ram
     Amount of RAM to request, in G (default = $ram).

   --cores
     CPU cores to request (default = $cores).

   --eddy-correct-method 
     Eddy correction method, either ANTs or FSL (default = $eddyCorrectMethod).
  

  Hard-coded settings:

  There are various hard-coded settings for software and templates. See the script and check these are correct.


  Output:

    See processScanDTI.pl for details of output. This script is a wrapper that sets up I/O and logs for running processScanDTI.pl.

};

if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

my $subject = "";
my @timepoints = ();

GetOptions ("subject=s" => \$subject,
	    "timepoints=s{1,1000}" => \@timepoints,
	    "dwi-base-dir=s" => \$dwiBaseDir,
	    "antsct-base-dir=s" => \$antsCTBaseDir,
	    "output-base-dir=s" => \$outputBaseDir,
	    "qsub=i" => \$submitToQueue,
	    "ram=s" => \$ram,
	    "cores=i" => \$cores,
            "eddy-correct-method=s" => \$eddyCorrectMethod
    )
    or die("Error in command line arguments\n");

my $qVmem = "${ram}G";

# set Camino memory limit
my $caminoHeapSize = $ram * 800;

my $numTimepoints = scalar(@timepoints);

if ($numTimepoints == 0) {
    @timepoints = `ls ${dwiBaseDir}/${subject}`;

    chomp @timepoints;

    $numTimepoints = scalar(@timepoints);

    if ($numTimepoints == 0) {
	print "\n  No DWI data to process for $subject\n";
	exit 1;
    }
}

my $parallel="";

if ($cores > 1) {
    $parallel="-pe unihost $cores -binding linear:$cores"; 
}

# Submit each time point

for ( my $i=0; $i < $numTimepoints; $i++ ) {
    my $tp=$timepoints[$i];
    
    my $inputDWI_Dir="${dwiBaseDir}/${subject}/${tp}/DWI";
    my $antsCTOutputRoot="${antsCTBaseDir}/${subject}/${tp}/${subject}_${tp}_";
    
    # Check for correct input
    if (! -f "${antsCTOutputRoot}CorticalThickness.nii.gz") {
	print "\n  Incomplete or missing ANTsCT data for $subject $tp \n";
	next;
    }

    if (! -d $inputDWI_Dir ) {
	print "\n  No DWI data for $subject $tp\n";
	next;
    }

    # Quick check for some data
    my @bvecs = `ls ${inputDWI_Dir}/*.bvec 2> /dev/null`;

    if (scalar(@bvecs) == 0) {
	print "\n  No DTI data for $subject $tp\n";
	next;
    }

    # Match acqParams to protocol used for this DWI data
    $bvecs[0] =~ m/${subject}_${tp}_[0-9]{4}_(.*).bvec/;

    my $dwiProtocol = $1;

    my $acqParams = "${acqParamsDir}/${dwiProtocol}_acqp.txt";
    
    my $outputTPDir="${outputBaseDir}/${subject}/${tp}";     
    
    if (! -d $outputTPDir) {
	mkpath($outputTPDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $outputTPDir\n\t";
    }
    if (! -d "${outputTPDir}/logs") { 	
      mkpath("${outputTPDir}/logs", {verbose => 0, mode => 0775}) or die "Cannot create output directory ${outputTPDir}/logs\n\t";
    }
    if (! -d "${outputTPDir}/scripts") {
      mkpath("${outputTPDir}/scripts", {verbose => 0, mode => 0775}) or die "Cannot create output directory ${outputTPDir}/scripts\n\t";
    }
    
    # increment log counter as needed
    my $runCounter=1;
    
    my $runCounterFormatted = sprintf("%03d", $runCounter);
    
    my $logFile="${outputTPDir}/logs/dti_${subject}_${tp}_log${runCounterFormatted}.txt";
    
    while (-f $logFile) {

	$runCounter = $runCounter + 1;
	
	$runCounterFormatted = sprintf("%03d", $runCounter);

	$logFile="${outputTPDir}/logs/dti_${subject}_${tp}_log${runCounterFormatted}.txt";
    }

    if ($runCounter > 1) {
        print "  Some output exists for ${subject} ${tp}, resubmitting\n";
    }

    my $scriptToRun="${outputTPDir}/scripts/dti_${subject}_${tp}_${runCounterFormatted}.sh";

    my $antsVersion = `cat ${antsPath}version.txt`;
    chomp $antsVersion;

    my $fh;

    open($fh, ">", $logFile);

    print $fh "ANTs version: ${antsVersion}\n\n";
    
    close($fh);

    open($fh, ">", $scriptToRun);
    
    # Set paths here so they can't get altered by user profiles
    
    # Just ANTs and Camino; FSL version is hard coded into eddy scripts
    
    print $fh qq{
export ANTSPATH=$antsPath

export PATH=${caminoDir}:${antsPath}:\${PATH}

export CAMINO_HEAP_SIZE=${caminoHeapSize};

${Bin}/processScanDTI.pl --input-dir $inputDWI_Dir --output-dir $outputTPDir --output-file-root ${subject}_${tp}_ --acq-params $acqParams --eddy-correct-method $eddyCorrectMethod --antsct-output-root $antsCTOutputRoot --template $template --standard-template $mniTemplate --standard-template-warp-root $templateToMNI_WarpRoot --template-mask $templateMask --standard-template-mask $mniTemplateMask
};

    close $fh;
    
    if ($submitToQueue) {
	system("qsub -l h_vmem=${qVmem},s_vmem=${qVmem} $parallel -cwd -S /bin/bash -j y -o $logFile $scriptToRun");
	system("sleep 0.25");
    }
    else {
	system("/bin/bash $scriptToRun >> $logFile 2>&1");
    }
    
    
} # for all timepoints



