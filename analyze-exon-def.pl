#!/usr/bin/env perl
use strict;
use ProgramName;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;
use TempFilename;
use Getopt::Std;
our $opt_q;
getopts('q');

my $slurm=$ENV{"SLURM_JOB_ID"};

my $name=ProgramName::get();
die "
$name [options] <config-file> <reference.multi-fasta> <local.gff> <max-VCF-errors>
  * local.gff must be a GTF/GFF2 file, not GFF3; all elements must
    have transcript_id and gene_id elements.  You can use the local.gff
    file produced by make-individual-genomes.pl
  * VCF files must be bgzipped (not gzipped) and indexed with tabix,
    may not contain uncalled variants (./.) or have extra attributes in
    the genotype field; column 7 must be PASS and column 9 must be GT
  * chromosome names must begin with the letters \"chr\" in all files

Please refer to http://geneprediction.org/ACE for detailed instructions.

" unless @ARGV==4;
my ($modelDir,$refFasta,$refGFF,$MAX_VCF_ERRORS)=@ARGV;
my $commandline = join " ", $0, @ARGV;

my $QUIET=$opt_q ? "-q" : "";
my $DEBUG=0; # THIS ALSO TERMINATES THE SCRIPT AFTER THE FIRST GENE
my $MAX_TRANSCRIPTS=-1;
my $STOP_AFTER; #="ENST00000250572.8";
my $ACE=$ENV{"ACEPLUS"};
if(length($ACE)<1) { system("env"); die "environment variable ACEPLUS is not set" }
my $refFastaTemp=TempFilename::generate();
my $oneGeneGFF=TempFilename::generate();
my $tempRevcomp=TempFilename::generate();
my $gffReader=new GffTranscriptReader();
my $byGene=$gffReader->loadGeneIdHash($refGFF);
my $refReader=new FastaReader($refFasta);
my $fastaWriter=new FastaWriter();

my $totalTranscripts=0;
while(1) {
  my ($refDef,$refSeqRef)=$refReader->nextSequenceRef();
  last unless $refDef;
  $refDef=~/^\s*>\s*(\S+)_(\d)/ 
    || die "Can't parse ID from defline: $refDef";
  my ($refID,$hap)=($1,$2);
  my $modelFile=getModelFile($refSeqRef,$modelDir);
  my $transcripts=$byGene->{$refID};
  next unless $transcripts && $transcripts->[0];
  my $strand=$transcripts->[0]->getStrand();
  $fastaWriter->writeFastaFromRef($refDef,$refSeqRef,$refFastaTemp);
  if($strand eq "-") {
    System("revcomp-fasta.pl $refFastaTemp > $tempRevcomp ; mv $tempRevcomp $refFastaTemp") }
  my $numTranscripts=@$transcripts;
  for(my $i=0 ; $i<$numTranscripts ; ++$i) {
    my $transcript=$transcripts->[$i];
    my $transID=$transcript->getTranscriptId();
    next unless $STOP_AFTER eq "" || $STOP_AFTER eq $transID;
    ++$totalTranscripts;
    if($MAX_TRANSCRIPTS>0 && $totalTranscripts>$MAX_TRANSCRIPTS) { exit }
    my $gff=$transcript->toGff();
    open(GFF,">$oneGeneGFF") || die "can't write file $oneGeneGFF";
    print GFF "$gff";
    close(GFF);
    if($strand eq "-") {
      my $substrateLength=length($$refSeqRef);
      System("revcomp-gff.pl $oneGeneGFF $substrateLength > $tempRevcomp ; mv $tempRevcomp $oneGeneGFF"); }
    my $errorsFlag=$MAX_VCF_ERRORS>=0 ? "-e $MAX_VCF_ERRORS" : "";
    my $command="$ACE/analyze-exon-def $errorsFlag $modelFile $oneGeneGFF $refFastaTemp";
    system($command);
    if($DEBUG || $STOP_AFTER eq $transID) {
      print "$command\n";
      print "exiting for debugging\n";
      exit;
    }
    last; # one transcript from each gene
  }
}
unlink($refFastaTemp);
unlink($oneGeneGFF); unlink($tempRevcomp);

sub System
{
  my ($cmd)=@_;
  print("$cmd\n") if $DEBUG;
  system($cmd);
}

sub getModelFile
{
  my ($seqRef,$dir)=@_;
  return $dir unless -d $dir;
  my $ATCG=$$seqRef=~s/([ACGT])/$1/g;
  my $GC=$$seqRef=~s/([CG])/$1/g;
  my $gc=$GC/$ATCG;
  my $rounded=int($gc*100+5/9)/100;
  my $range;
  if($gc<=0.43) { $range="0-43" }
  elsif($gc<=0.51) { $range="43-51" }
  elsif($gc<=0.57) { $range="51-57" }
  else { $range="57-100" }
  return "$dir/ace.$range.config";
}






