
THIS FILE IS OBSOLETE


#!/bin/env perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <infile.gff> <outfile.gff> <margin>" unless @ARGV==3;
my ($infile,$outfile,$MARGIN)=@ARGV;

open(OUT,">$outfile") || die $outfile;
my $reader=new GffReader;
while(1) {
  my $feature=$reader->nextFeature();
  last unless $feature;
  my $begin=$feature->getBegin();
  my $numTrans=$gene->getNumTranscripts();
  for(my $i=0 ; $i<$numTrans ; ++$i) {
    my $transcript=$gene->getIthTranscript($i);
    $feature->shiftCoords($MARGIN-$begin);
  }
  my $gff=$feature->toGff();
  print OUT $gff;
}
close(OUT);

