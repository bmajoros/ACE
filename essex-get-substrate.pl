#!/usr/bin/env perl
use strict;
use EssexParser;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.essex> <substrateID>\n" unless @ARGV==2;
my ($infile,$ID)=@ARGV;

my $parser=new EssexParser($infile);
while(1) {
  my $report=$parser->nextElem();  last unless $report;
  my $thisID=$report->getAttribute("substrate");
  if($thisID eq $ID) {
    $report->print(\*STDOUT);
    print "\n";
  }
}


