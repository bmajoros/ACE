#!/bin/env perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <gencode.gff> <out.gff>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my @keepFeatures=("exon","UTR","CDS","transcript");
my %keepFeatures; foreach my $t (@keepFeatures) { $keepFeatures{$t}=1 }
my @geneTypes=("protein_coding","lincRNA");
my %geneTypes; foreach my $t (@geneTypes) { $geneTypes{$t}=1 }

my $reader=new GffReader();
open(IN,$infile) || die "can't open file $infile\n";
open(OUT,">$outfile") || die "can't write to file $outfile\n";
while(1) {
  my $feature=$reader->nextFeature(\*IN);
  last unless $feature;
  my $featureType=$feature->getType();
  next unless $keepFeatures{$featureType};
  my $extra=$feature->getExtraAsHash();
  my $transcriptType=$extra->{"transcript_type"};
  next unless $geneTypes{$transcriptType};
  next if $transcriptType eq "protein_coding" && $featureType eq "exon";
  my $status=$extra->{"transcript_status"};
  next unless $status eq "KNOWN";
  my $geneID=$extra->{"gene_id"};
  my $transID=$extra->{"transcript_id"};
  my $substrate=$feature->getSubstrate();
  my $source=$feature->getSource();
  my $begin=$feature->getBegin()+1; # FIXED ON 4/21/2016
  my $end=$feature->getEnd();
  my $score=$feature->getScore();
  my $strand=$feature->getStrand();
  my $frame=$feature->getFrame();
  print OUT "$substrate\t$source\t$featureType\t$begin\t$end\t$score\t$strand\t$frame\tgene_id=$geneID;transcript_id=$transID;type=$transcriptType;\n";
}
close(OUT);
close(IN);




