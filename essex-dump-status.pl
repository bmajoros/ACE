#!/bin/env perl
use strict;
use EssexParser;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.essex>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $parser=new EssexParser($infile);
while(1) {
  my $report=$parser->nextElem();
  last unless $report;
  my $transcriptID=$report->getAttribute("transcript-ID");
  my $substrate=$report->getAttribute("substrate");
  my $warnings=$report->getAttribute("vcf-warnings");
  my $errors=$report->getAttribute("vcf-errors");
  my $status=$report->findChild("status");
  next unless $status;
  my $code=$status->getIthElem(0);
  next unless $code;
  print "$substrate\t$warnings/$errors\t$transcriptID\t$code";
  my $numElem=$status->numElements();
  for(my $i=1 ; $i<$numElem ; ++$i) {
    my $elem=$status->getIthElem($i);
    if(EssexNode::isaNode($elem)) {
      my $tag=$elem->getTag();
      print " $tag";
    }
    else { print " $elem" }
  }
  print "\n";
}


