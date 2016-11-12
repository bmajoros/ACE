#!/bin/env perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <cigar-string>\n" unless @ARGV==1;
my ($cigar)=@ARGV;

print "ref\talt\tevent\n";
my $ref=0; my $alt=0;
while(length($cigar)>0) {
  $cigar=~/^(\d+)([MID])(.*)/ || die $cigar;
  my ($L,$op,$rest)=($1,$2,$3);
  if($op eq "M") { $ref+=$L; $alt+=$L }
  elsif($op eq "I") { print "$ref\t$alt\t${L} bp insertion\n"; $alt+=$L }
  elsif($op eq "D") { print "$ref\t$alt\t${L} bp deletion\n"; $ref+=$L }
  $cigar=$rest;
}


