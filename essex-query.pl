#!/usr/bin/env perl
use strict;
use EssexParser;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <in.essex> <query>\n" unless @ARGV==2;
my ($infile,$query)=@ARGV;

my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem();
  last unless $root;
  if(EssexNode::isaNode($root)) {
    if($root->getTag() eq $query) {
      $root->print(\*STDOUT);
      print "\n";
      next
    }
    my $array=$root->findDescendents($query);
    foreach my $elem (@$array) {
      $elem->print(\*STDOUT);
      print "\n";
    }
  }
}
