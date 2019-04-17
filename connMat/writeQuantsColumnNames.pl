#!/usr/bin/perl -w


use strict;

# Simple script to take in a list of nodes or whatever and write out a QuANTs compatible label definition

my $corticalLabelFile = "mindBoggleCorticalGraphNodes.csv";

my $fh;

open($fh, "<", $corticalLabelFile);

my $oldHeader = <$fh>;

my $newHeader = "Label.ID,QuantsLabelName,LongLabelName";

print "$newHeader\n";

while (my $line = <$fh>) {
  chomp($line);	

  my @tokens = split(",", $line);

  my $intensity = $tokens[0];

  print "${intensity},mindboggle_${intensity},$tokens[1]\n";
  
}

close($fh);
