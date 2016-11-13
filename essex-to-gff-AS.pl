#!/usr/bin/env perl
use strict;
use EssexParser;
use EssexACE;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <in.essex> <out.gff> <allele#>\n" unless @ARGV==3;
my ($infile,$outfile,$hap)=@ARGV;

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
  my $status=$ace->getStatusString();
  #if($status->hasDescendentOrDatum("bad-annotation")) { next }
  #if($status eq "mapped") {
  my $transcript=$ace->getMappedTranscript();
  if($transcript) {
    my $id=$transcript->getTranscriptId();
    $id.="_$hap";
    $transcript->setTranscriptId($id);
    $id=$transcript->getGeneId();
    $id.="_$hap";
    $transcript->setGeneId($id);
    print OUT $transcript->toGff();
  }
  if($status eq "splicing-changes") {
    my $transcripts=$ace->getAltTranscripts();
    my $n=@$transcripts;
    for(my $i=0 ; $i<$n ; ++$i) {
      my $transcript=$transcripts->[$i];
      my $id=$transcript->getTranscriptId();
      if($transcript->getSource() eq "SIMULATION") { $id="SIM$i\_$id\_$hap" }
      else { $id="ALT$i\_$id\_$hap" }
      $transcript->setTranscriptId($id);
      $id=$transcript->getGeneId();
      $id.="_$hap";
      $transcript->setGeneId($id);
      print OUT $transcript->toGff();
    }
  }
  elsif($status eq "mapped") {
    my $transcriptNode=
      $root->pathQuery("report/status/new-upstream-start-codon/transcript");
    if($transcriptNode) {
      my $transcript=new Transcript($transcriptNode);
      my $id=$transcript->getTranscriptId();
      $id="NEWSTART_$id\_$hap";
      $transcript->setTranscriptId($id);
      $id=$transcript->getGeneId();
      $id.="_$hap";
      $transcript->setGeneId($id);
      print OUT $transcript->toGff();
    }
    undef $transcriptNode;
  }
  undef $root; undef $ace;
}
close(OUT);
print "[done]\n";

