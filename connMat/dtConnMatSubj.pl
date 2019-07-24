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

my $antsPath = "/data/grossman/pcook/lausanneConMat/bin/ants/";

# Path where ANTsR etc can be found
my $rLibsUserPath = "/data/grossman/pipedream2018/bin/R/R-3.4.3/library";

my $caminoDir = "/data/grossman/pcook/lausanneConMat/bin/camino/bin";

# Two cores because we have Java and R running at the same time to reduce disk I/O
my $cores=2;

my $ram="8";
my $submitToQueue=1;

# CSV file containing JLF cortical label definitions
my $jlfCorticalLabelDef = "${Bin}/mindBoggleCorticalGraphNodes.csv";

# These are what are used to label the cortical mask from JLF, and become the graph nodes
# Default to the JLF labels themselves
my $graphCorticalLabelSystem = "PG_antsLabelFusionLabels";
my $graphCorticalLabelDef = $jlfCorticalLabelDef;


my $usage = qq{

  $0  
      --subject
      --antsct-base-dir
      --dt-base-dir
      [ options ]

  Required args:

   --subject
     Subject ID

   --antsct-base-dir | --antslongct-base-dir  
     Base ants[long]CT dir for T1 data. If longitudinal, set --longitudinal-target.

   --dt-base-dir
     Base DTI dir. There should be DT data for the time point(s) to be processed. For longitudinal data,
     there should also be FA images in the single-subject template space.

   --output-base-dir
     Base output directory. Output is organized under this as subj/tp/connMat_graphCorticalLabelSystem.

 Options:

   --timepoints  
     Timepoint(s) to process. Default is to process all available timepoints for the subject.

   --longitudinal-target
     Either "session" to evaluate connectivity in the intra-session T1, or "sst", to do it in the SST space.  
     If "session", labels for node creation must be defined for each time point.
     Default is "sst" if "--antslongct-base-dir" is present.

   --cortical-label-system
     System for labeling cortex. There must be an image in the antsCT output directory named subj_tp_system.nii.gz,
     containing cortical labels. By default, this is the JLF labels. 

     If this option is specified, the list of cortical label IDs and names must be provided with the 
     --cortical-label-def option. The image may contain non-cortical labels, but only the cortical labels will be 
     used for the graph (default = $graphCorticalLabelSystem).

   --cortical-label-def
     A CSV file containing "Label.ID,Label.Name", with cortical labels only (default = $graphCorticalLabelDef).

   --qsub
     Submit processing jobs to qsub (default = $submitToQueue).

   --ram
     Amount of RAM to request, in G (default = $ram). Increase this for for sub-1mm T1 scans.

   --cores 
     Number of CPU cores to request (default = $cores).


  Wrapper script for building matrices of connectivity. This script just submits jobs to the queue.

  The target T1 image / SST should be labeled using the pseudo-geodesic JLF.

  DT data is read from

    \${dt-base-dir}/subj/tp/

  The DTI preprocessing pipeline should have been run, such that dt/ and distCorr/ exist; 
  these will be used to do the tracking and transfer the results to T1 space.

  
  Output:

    See dtConnMatCortical.pl and dtConnMat.pl.

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
my @timepoints = ();

# If empty string, expect cross-sectional antsct input
my $longitudinalTarget = "";


GetOptions ("subject=s" => \$subject,
	    "timepoint|timepoints=s{1,1000}" => \@timepoints,
	    "antsct-base-dir=s" => \$antsCTBaseDir,
	    "antslongct-base-dir=s" => sub {
		my ($opt_name, $opt_value) = @_; 
		$antsCTBaseDir = $opt_value; 
		if ($longitudinalTarget eq "") {
		    $longitudinalTarget = "sst";
		}
	    },
	    "dt-base-dir=s" => \$dtBaseDir,
	    "longitudinal-target=s" => sub { 
		my ($opt_name, $opt_value) = @_; 
		$longitudinalTarget = lc($opt_value); 
	    },
	    "output-base-dir=s" => \$outputBaseDir,
	    "cortical-label-system=s" => \$graphCorticalLabelSystem,
	    "cortical-label-def=s" => \$graphCorticalLabelDef,
	    "qsub=i" => \$submitToQueue,
	    "ram=s" => \$ram,
	    "cores=i" => \$cores
    )
    or die("Error in command line arguments\n");

my $numTimepoints = scalar(@timepoints);

if ($numTimepoints == 0) {
    @timepoints = `ls ${dtBaseDir}/${subject} | grep -v singleSubjectTemplate`;

    chomp @timepoints;

    $numTimepoints = scalar(@timepoints);

    if ($numTimepoints == 0) {
	print "\n  No DT data to process for $subject\n";
	exit 1;
    }
}

if (! -d $outputBaseDir) {
    mkpath($outputBaseDir, {verbose => 0, mode => 0775}) or die "Output base directory \n  $outputBaseDir \n does not exist and cannot be created \n\t";
}
 
my $qVmem = "${ram}G";

# set Camino memory limit
my $caminoHeapSize = 3200;

my $parallel="";

if ($cores > 1) {
    $parallel="-pe unihost $cores -binding linear:$cores"; 
}


# Submit each time point

for ( my $i=0; $i < $numTimepoints; $i++ ) {
    my $tp = $timepoints[$i];

    # Base dir containing DT stuff for this TP
    my $tpDTIDir = "${dtBaseDir}/${subject}/${tp}";

    my $dt = "${tpDTIDir}/dt/${subject}_${tp}_DT.nii.gz";

    if (! -f "$dt" ) {
	print "\n  No DT for $subject $tp\n";
	next;
    }

    my $tpOutputDir = "${outputBaseDir}/${subject}/${tp}/connMat_${graphCorticalLabelSystem}";     

    if (-d $tpOutputDir) {
	print "\n  Output dir for system ${graphCorticalLabelSystem} already exists for $subject $tp \n";
	next;
    }

    mkpath($tpOutputDir, {verbose => 0, mode => 0775}) or die "Cannot create output directory $tpOutputDir\n\t";
    
    my $antsCTArgString = "";

    if ($longitudinalTarget ne "") {
	$antsCTArgString = "--antslongct-base-dir ${antsCTBaseDir} \\\n    --longitudinal-target $longitudinalTarget";
    }
    else {
	$antsCTArgString = "--antsct-base-dir $antsCTBaseDir";
    }
    

    # Root for output we will create
    my $tpOutputRoot = "${tpOutputDir}/${subject}_${tp}_";

    my $logFile = "${tpOutputDir}/connMat_${subject}_${tp}_log.txt";
    
    my $scriptToRun = "${tpOutputDir}/connMat_${subject}_${tp}.sh";

    my $fh;

    open($fh, ">", $scriptToRun);
    
    # Set paths here so they can't get altered by user profiles
    print $fh qq{
export ANTSPATH=${antsPath}

export PATH=${caminoDir}:${antsPath}:\${PATH}

export CAMINO_HEAP_SIZE=${caminoHeapSize}

export R_LIBS=""

export R_LIBS_USER=${rLibsUserPath}

${Bin}/dtConnMatCortical.pl \\
    --subject $subject \\
    --timepoint $tp \\
    --dt-base-dir $dtBaseDir \\
    $antsCTArgString \\
    --cortical-label-system $graphCorticalLabelSystem \\
    --cortical-label-def $graphCorticalLabelDef \\
    --output-base-dir $outputBaseDir 

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

