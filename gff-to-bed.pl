#!/usr/bin/env perl
use strict;
use GffTranscriptReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff> <out.bed>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write file: $outfile\n";
my $reader=new GffTranscriptReader();
my $genes=$reader->loadGenes($infile);
my $n=@$genes;
for(my $i=0 ; $i<$n ; ++$i) {
  my $gene=$genes->[$i];
  my $chr=$gene->getSubstrate();
  my $begin=$gene->getBegin();
  my $end=$gene->getEnd();
  my $strand=$gene->getStrand();
  my $id=$gene->getId();
  print OUT "$chr\t$begin\t$end\t$id\t0\t$strand\n";
}
close(OUT);



