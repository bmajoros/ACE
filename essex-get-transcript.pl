#!/usr/bin/env perl
use strict;
use EssexParser;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.essex> <transcriptID>\n" unless @ARGV==2;
my ($infile,$ID)=@ARGV;

my $parser=new EssexParser($infile);
while(1) {
  my $report=$parser->nextElem();  last unless $report;
  my $transcriptID=$report->getAttribute("transcript-ID");
  if($transcriptID eq $ID) {
    $report->print(\*STDOUT);
    print "\n";
    last;
  }
}


