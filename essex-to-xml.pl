#!/bin/env perl
use strict;
use ProgramName;
use EssexParser;

my $name=ProgramName::get();
die "$name <in.sx>\n" unless @ARGV>=1;
my ($infile)=@ARGV;

print "<ace>\n";
my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem(); last unless $root;
  $root->printRecursiveXML(1,\*STDOUT);
  print "\n";
}
print "</ace>\n";


