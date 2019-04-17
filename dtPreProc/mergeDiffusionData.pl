#!/usr/bin/perl -w
#
# Merge DWI data and bvecs for further processing
#


use strict;
use FindBin qw($Bin);
use File::Copy;
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
     --input-dir
     --output-root

  Required args:

   --input-dir
     Directory to search for diffusion data, in the format root.bvec, root.bval, root.nii.gz. 
     Other files without a matching .bvec are ignored.
  
   --output-root
     Output root, prepended to output file names.

   Requires ANTSPATH. Current ANTSPATH is : $antsPath

   Produces outputRoot[.nii.gz,.bval,.bvec] containing concatenated DWI data.

   Series are concatenated in the order they are listed with ls. The ordering determines the choice of reference 
   space for eddy correction and for the distortion correction to T1. 

   There must be bvecs, bvals, and image data for each series. For example, name.bvec, name.bval and name.nii.gz.

};


if ($#ARGV < 0) {
    print $usage;
    exit 1;
}

my ($inputDir, $outputRoot) = @ARGV;;



GetOptions ("input-dir=s" => \$inputDir,   
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


# Search for DWI data, in ${inputDir} only
my @bvecs = `ls -d $inputDir/* | grep ".bvec"`; 

chomp @bvecs;

my $numScans = scalar(@bvecs);

my @bvals = ();
my @dwi = ();

for (my $i = 0; $i < $numScans; $i++) {
    my $image = $bvecs[$i];
    
    $image =~ s/\.bvec/.nii.gz/;

    if (! -f $image) {
	die("\n  No image matching ${bvecs[$i]} ");
    }

    my $bval = $bvecs[$i];

    $bval =~ s/\.bvec/.bval/;

    if (! -f $bval) {
	die("\n  No bval file matching ${bvecs[$i]} ");
    }
    
    push(@bvals, $bval);
    push(@dwi, $image);
}

if ($numScans == 0) {
    print "\n  No DWI data found in $inputDir \n";
    exit 1;
}

# Directory for temporary files that is deleted later if $cleanup
my $tmpDir = "";

my $tmpDirBaseName = "${outputFileRoot}mergeDWI";

if ( !($sysTmpDir && -d $sysTmpDir) ) {
    $tmpDir = "${outputDir}/${tmpDirBaseName}";
}
else {
    # Have system tmp dir
    $tmpDir = "${sysTmpDir}/${tmpDirBaseName}";
}

# Gets removed later, so check we can create this and if not, exit immediately
mkpath($tmpDir, {verbose => 0, mode => 0755}) or die "Cannot create working directory $tmpDir (maybe it exists from a previous failed run)\n\t";


print "\nProcessing " . $numScans . " scans.\n";

print "\nScan merge order\n:";

foreach my $file (@bvals) {

    my $baseName = basename($file, (".bval"));
    
    print "\t${baseName}\n";
} 


my $bvalMaster = "${tmpDir}/${outputFileRoot}.bval";
my $bvecMaster = "${tmpDir}/${outputFileRoot}.bvec";

if ( $numScans > 1 ) {
    my $bvalFileNames = join(" ", @bvals);

    system("paste -d \" \" $bvalFileNames > $bvalMaster");

    my $fh;

    # Make sure there is no double spacing - it breaks space delimited ITK readers
    open($fh, "<", "$bvalMaster") or die "Cant open $bvalMaster\n";

    local $/ = undef;

    my $bvalString = <$fh>;

    close($fh);

    $bvalString =~ s/[ ]{2,}/ /g;

    open($fh, ">", "$bvalMaster") or die "Cant open $bvalMaster\n";

    print $fh $bvalString;

    close($fh);

    my $bvecFileNames = join(" ", @bvecs);

    system("paste -d \" \" $bvecFileNames > $bvecMaster");

    # Make sure there is no double spacing - it breaks space delimited ITK readers
    open($fh, "<", "$bvecMaster") or die "Cant open $bvecMaster\n";

    local $/ = undef;

    my $bvecString = <$fh>;

    close($fh);

    $bvecString =~ s/[ ]{2,}/ /g;

    open($fh, ">", "$bvecMaster") or die "Cant open $bvecMaster\n";

    print $fh $bvecString;

    close($fh);
}
else {
    copy($bvals[0], $bvalMaster);
    copy($bvecs[0], $bvecMaster);
}

my $mergedDWI = "${tmpDir}/${outputFileRoot}.nii.gz";

copy($dwi[0], $mergedDWI);

if ( $numScans > 1 ) {
    for ( my $i = 1; $i < scalar(@dwi); $i += 1) {
        system("${antsPath}/ImageMath 4 $mergedDWI stack $mergedDWI $dwi[$i]");
    }
}

# Check merge happened correctly
my $numMeasExpected = `wc -w $bvalMaster | cut -d ' ' -f 1`;

my $dwiMasterInfo = `${antsPath}/PrintHeader $mergedDWI`;

$dwiMasterInfo =~ m/dim\[4\] = (\d+)/;

my $masterComps = $1;

if ($masterComps == 1) {
  # check 5th dimension
  $dwiMasterInfo =~ m/dim\[5\] = (\d+)/;
  $masterComps = $1;
}

if ($masterComps != $numMeasExpected) {
  die ("\nDWI data is inconsistent with number of bvalues, expected $numMeasExpected volumes but found $masterComps. Either the data could not be merged or the bvals are incorrect\n");
}

system("mv $mergedDWI $outputDir");
system("mv $bvalMaster $outputDir");
system("mv $bvecMaster $outputDir");


if ($cleanup) {
    system("rm -f ${tmpDir}/*");
    system("rmdir ${tmpDir}");
}
