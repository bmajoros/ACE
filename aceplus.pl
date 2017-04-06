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
if($slurm) { print "SLURM $slurm\n" }

my $name=ProgramName::get();
die "
$name [options] <config-file> <reference.multi-fasta> <personal.multi-fasta> <local.gff> <max-VCF-errors> <out.essex>

  -q = quiet: don't report annotation errors or transcripts that map perfectly
  * local.gff must be a GTF/GFF2 file, not GFF3; all elements must
    have transcript_id and gene_id elements.  You can use the local.gff
    file produced by make-individual-genomes.pl
  * VCF files must be bgzipped (not gzipped) and indexed with tabix,
    may not contain uncalled variants (./.) or have extra attributes in
    the genotype field; column 7 must be PASS and column 9 must be GT
  * chromosome names must begin with the letters \"chr\" in all files

Please refer to http://geneprediction.org/ACE for detailed instructions.

" unless @ARGV==6;
my ($modelDir,$refFasta,$altFasta,$refGFF,$MAX_VCF_ERRORS,$outACE)=@ARGV;
my $commandline = join " ", $0, @ARGV;
print "$commandline\n";
system("hostname");

my $QUIET=$opt_q ? "-q" : "";
my $DEBUG=0; # THIS ALSO TERMINATES THE SCRIPT AFTER THE FIRST GENE
my $MAX_TRANSCRIPTS=-1;
my $STOP_AFTER; #="ENST00000250572.8";
my $START_AT="ENSG00000231982.1";
my $ACE=$ENV{"ACEPLUS"};
if(length($ACE)<1) { system("env"); die "environment variable ACEPLUS is not set" }
my $refFastaTemp=TempFilename::generate();
my $altFastaTemp=TempFilename::generate();
#my $labelingTemp=TempFilename::generate();
my $oneGeneGFF=TempFilename::generate();
my $aceTemp=TempFilename::generate();
my $gffTemp=TempFilename::generate();
#my $tempXML=TempFilename::generate();
my $tempRevcomp=TempFilename::generate();
my $gffReader=new GffTranscriptReader();
my $byGene=$gffReader->loadGeneIdHash($refGFF);
my $refReader=new FastaReader($refFasta);
my $altReader=new FastaReader($altFasta);
my $fastaWriter=new FastaWriter();
#unlink($outLab) if -e $outLab;
unlink($outACE) if -e $outACE;
#unlink($outXML) if -e $outXML;

system("date");
my $totalTranscripts=0;
#open(LABELING,">$outLab") || die "Can't write to $outLab";
while(1) {
  my ($altDef,$altSeqRef)=$altReader->nextSequenceRef();
  last unless $altDef;
  $altDef=~/^\s*>\s*(\S+)_(\d)/ 
    || die "Can't parse ID from alt defline: $altDef";
  my ($altID,$hap)=($1,$2);
  if($START_AT ne "") {
    if($altID eq $START_AT) { $START_AT="" }
    else { next }
  }
  print "$altID\_$hap\n";
  my $modelFile=getModelFile($altSeqRef,$modelDir);
  my $transcripts=$byGene->{$altID};
  #die "no transcripts on substrate $altID" 
  next
    unless $transcripts && $transcripts->[0];
  my $strand=$transcripts->[0]->getStrand();
  my ($refDef,$refSeqRef);
  while(1) {
    ($refDef,$refSeqRef)=$refReader->nextSequenceRef();
    die "no more sequences in $refFasta\n" unless $refDef;
    $refDef=~/^\s*>\s*(\S+)/ || die "Can't parse ID from ref defline: $refDef";
    my $refID=$1;
    if($refID=~/(\S+)_\d+$/) { $refID=$1 }
    last if $refID eq $altID;
  }
  $fastaWriter->writeFastaFromRef($refDef,$refSeqRef,$refFastaTemp);
  $fastaWriter->writeFastaFromRef($altDef,$altSeqRef,$altFastaTemp);
  if($strand eq "-") {
    System("revcomp-fasta.pl $refFastaTemp > $tempRevcomp ; mv $tempRevcomp $refFastaTemp");
    System("revcomp-fasta.pl $altFastaTemp > $tempRevcomp ; mv $tempRevcomp $altFastaTemp");
  }
  my $numTranscripts=@$transcripts;
  for(my $i=0 ; $i<$numTranscripts ; ++$i) {
    my $transcript=$transcripts->[$i];
    my $transID=$transcript->getTranscriptId();
    print "transcript $transID\n";
    next unless $STOP_AFTER eq "" || $STOP_AFTER eq $transID;
    ++$totalTranscripts;
    if($MAX_TRANSCRIPTS>0 && $totalTranscripts>$MAX_TRANSCRIPTS) { exit }
    my $gff=$transcript->toGff();
    open(GFF,">$oneGeneGFF") || die "can't write file $oneGeneGFF";
    print GFF "$gff";
    close(GFF);
    my $revFlag="";
    if($strand eq "-") {
      my $substrateLength=length($$refSeqRef);
      System("revcomp-gff.pl $oneGeneGFF $substrateLength > $tempRevcomp ; mv $tempRevcomp $oneGeneGFF");
      $revFlag="-c";
    }
    my $errorsFlag=$MAX_VCF_ERRORS>=0 ? "-e $MAX_VCF_ERRORS" : "";
    my $command="$ACE/aceplus $QUIET $errorsFlag $revFlag $modelFile $oneGeneGFF $refFastaTemp $altFastaTemp $gffTemp $aceTemp";
    print "$command\n" ;#if $DEBUG;
    my $err=`$command`;
    print "$err\n";
    if($err=~/abort/ || $err=~/INTERNAL ERROR/ || $err=~/no transcripts/
       || $err=~/no such file/) { exit }
    die "ACE terminated abnormally:\n$command" 
      unless $err=~/ACE terminated successfully/;
    System("cat $aceTemp >> $outACE");
    #System("cat $tempXML >> $outXML");
    #my $labeling;
    #open(IN,$labelingTemp) || die "Can't read $labelingTemp";
    #while(<IN>) { chomp; $labeling.=$_ }
    #close(IN);
    #print LABELING "$altID\n$labeling\n";
    if($DEBUG || $STOP_AFTER eq $transID) {
      print "$command\n";
      print "exiting for debugging\n";
      exit;
    }
  }
}
unlink($refFastaTemp); unlink($altFastaTemp); 
#unlink($labelingTemp);
unlink($aceTemp); unlink($oneGeneGFF); unlink($tempRevcomp); unlink($gffTemp);
#unlink($tempXML);
#close(LABELING);
close(ACE);
system("date");
print "[done]\n";

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






