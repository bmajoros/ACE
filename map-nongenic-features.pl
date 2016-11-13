#!/usr/bin/env perl
use strict;
use ProgramName;
use FastaReader;
use GffReader;
use TempFilename;

# Parse the command line
my $name=ProgramName::get();
die "$name <features.gff> <regions.fasta> <out.gff>\n" unless @ARGV==3;
my ($refGFF,$altFasta,$outfile)=@ARGV;

# Declare some globals
my $DEBUG=1;
my $ICE=$ENV{"ICE"};
my $cigarFile=TempFilename::generate();
my $gffTemp=TempFilename::generate();
my $mappedTemp=TempFilename::generate();
my $fastaReader=new FastaReader($altFasta);
my $gffReader=new GffReader();

# Load all features and hash by substrate
my $bySubstrate=$gffReader->hashBySubstrate($refGFF);

# Go sequentially through the FASTA file
while(1) {
  my ($altDef,$altSeqRef)=$fastaReader->nextSequenceRef();
  last unless $altDef;
  $altDef=~/^\s*>\s*(\S+)_(\d)/ 
    || die "Can't parse ID from alt defline: $altDef";
  my ($altID,$hap)=($1,$2);
  my $features=$bySubstrate->{$altID};
  next unless $features && @$features>0;
  $altDef=~/\/cigar=(\S+)/ || die "can't parse cigar string from defline: $altDef\n";
  my $cigar=$1;
  open(CIGAR,">$cigarFile") || die "Can't write file $cigarFile\n";
  print CIGAR "$cigar\n";
  close(CIGAR);
  open(TEMP,">gffTemp") || die "Can't write file $gffTemp\n";
  my $numFeatures=@$features;
  for(my $i=0 ; $i<$numFeatures ; ++$i) {
    my $feature=$features->[$i];
    my $gff=$feature->toGff();
    print TEMP $gff;
    my $command="$ICE/map-annotations $gffTemp $cigarFile $mappedTemp";
    System($command);
    System("cat $mappedTemp >> $outfile");
  }
}
unlink($cigarFile);
unlink($gffTemp);
unlink($mappedTemp);


sub System
{
  my ($cmd)=@_;
  print("$cmd\n") if $DEBUG;
  system($cmd);
}






