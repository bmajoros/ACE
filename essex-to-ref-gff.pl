#!/usr/bin/env perl
use strict;
use EssexParser;
use EssexACE;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <in.essex> <out.gff>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my %seen;
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem();
  last unless $root;
  my $ace=new EssexACE($root);
  my $transcriptID=$ace->getTranscriptID();
  next if $seen{$transcriptID};
  $seen{$transcriptID}=1;
  my $substrate=$root->getAttribute("substrate");
  #my $status=$ace->getStatusString();
  #if($status->hasDescendentOrDatum("bad-annotation")) { next }
  #next unless if($status eq "mapped");
  my $transcript=$ace->getRefTranscript();
  if($transcript) {
    my $id=$transcript->getTranscriptId();
    $transcript->setTranscriptId($id);
    $transcript->setSubstrate($substrate);
    $id=$transcript->getGeneId();
    $transcript->setGeneId($id);
    print OUT $transcript->toGff();
  }
  undef $root; undef $ace;
}
close(OUT);


