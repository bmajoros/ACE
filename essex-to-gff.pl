#!/usr/bin/env perl
use strict;
use EssexParser;
use EssexACE;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <in.essex> <out.gff>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem();
  last unless $root;
  my $ace=new EssexACE($root);
  my $transcriptID=$ace->getTranscriptID();
  my $status=$ace->getStatusString();
  if($status eq "mapped") {
    my $transcript=$ace->getMappedTranscript();
    print OUT $transcript->toGff();
  }
  elsif($status eq "splicing-changes") {
    my $transcripts=$ace->getAltTranscripts();
    my $n=@$transcripts;
    for(my $i=0 ; $i<$n ; ++$i) {
      my $transcript=$transcripts->[$i];
      my $id=$transcript->getTranscriptId();
      $id="ALT$i\_$id";
      $transcript->setTranscriptId($id);
      print OUT $transcript->toGff();
    }
  }
  undef $root; undef $ace;
}
close(OUT);


