#!/bin/env perl
use strict;
use EssexParser;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.essex>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my (%counts,%flags);
my $parser=new EssexParser($infile);
while(1) {
  my $report=$parser->nextElem();
  last unless $report;
  my $status=$report->findChild("status");
  next unless $status;
  my $code=$status->getIthElem(0);
  next unless $code;
  ++$counts{$code};
  my $numElem=$status->numElements();
  for(my $i=1 ; $i<$numElem ; ++$i) {
    my $elem=$status->getIthElem($i);
    if(EssexNode::isaNode($elem)) {
      my $tag=$elem->getTag();
      ++$flags{$code}->{$tag};
    }
    else { ++$flags{$code}->{$elem} }
  }
}

my @keys=keys %counts;
foreach my $key (@keys) {
  my $count=$counts{$key};
  print "$key\t$count\n";
  my @keys=keys %{$flags{$key}};
  foreach my $key2 (@keys) {
    my $count=$flags{$key}->{$key2};
    print "\t$key2\t$count\n";
  }
}


