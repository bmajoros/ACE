#!/usr/bin/env perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <chrom-name> <chrom-len> <bin-size> <out.bed>\n" unless @ARGV==4;
my ($chr,$chromLen,$binSize,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my $begin=0;
while($begin<$chromLen) {
  my $end=$begin+$binSize;
  if($end>$chromLen) { $end=$chromLen }
  print OUT "$chr\t$begin\t$end\n";
  $begin=$end;
}
close(OUT);



