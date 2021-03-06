#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;
use TempFilename;

my $slurm=$ENV{"SLURM_JOB_ID"};
print "SLURM $slurm\n";

my $name=ProgramName::get();
die "$name <model-file-or-dir> <ref.multi-fasta> <alt.multi-fasta> <ref.multi-gff> <out.ace>\n" unless @ARGV==5;
my ($modelDir,$refFasta,$altFasta,$refGFF,$outACE)=@ARGV;
my $commandline = join " ", $0, @ARGV;
print "$commandline\n";
system("hostname");

my $QUIET=""; #"-q";
my $DEBUG=0;
my $MAX_TRANSCRIPTS=-1;
my $MAX_VCF_ERRORS=0;
my $STOP_AFTER;#="ENST00000526104";
my $ACE=$ENV{"ACE"};
print "TRACE ", __LINE__, "\n";
my $refFastaTemp=TempFilename::generate();
my $altFastaTemp=TempFilename::generate();
#my $labelingTemp=TempFilename::generate();
my $oneGeneGFF=TempFilename::generate();
my $aceTemp=TempFilename::generate();
my $gffTemp=TempFilename::generate();
#my $tempXML=TempFilename::generate();
my $tempRevcomp=TempFilename::generate();
my $gffReader=new GffTranscriptReader();
print "TRACE ", __LINE__, "\n";
my $byGene=$gffReader->loadGeneIdHash($refGFF);
print "TRACE ", __LINE__, "\n";
my $refReader=new FastaReader($refFasta);
my $altReader=new FastaReader($altFasta);
my $fastaWriter=new FastaWriter();
#unlink($outLab) if -e $outLab;
print "TRACE ", __LINE__, "\n";
unlink($outACE) if -e $outACE;
#unlink($outXML) if -e $outXML;

system("date");
my $totalTranscripts=0;
#open(LABELING,">$outLab") || die "Can't write to $outLab";
while(1) {
print "TRACE ", __LINE__, "\n";
  my ($altDef,$altSeqRef)=$altReader->nextSequenceRef();
  last unless $altDef;
  $altDef=~/^\s*>\s*(\S+)_(\d)/ 
    || die "Can't parse ID from alt defline: $altDef";
  my ($altID,$hap)=($1,$2);
  print "$altID\_$hap\n";
  my $modelFile=getModelFile($altSeqRef,$modelDir);
print "TRACE ", __LINE__, "\n";
  my $transcripts=$byGene->{$altID};
  #die "no transcripts on substrate $altID" 
  next
    unless $transcripts && $transcripts->[0];
  my $strand=$transcripts->[0]->getStrand();
  my ($refDef,$refSeqRef);
  while(1) {
print "TRACE ", __LINE__, "\n";
    ($refDef,$refSeqRef)=$refReader->nextSequenceRef();
print "TRACE ", __LINE__, "\n";
    die "no more sequences in $refFasta\n" unless $refDef;
    $refDef=~/^\s*>\s*(\S+)/ || die "Can't parse ID from ref defline: $refDef";
    my $refID=$1;
    if($refID=~/(\S+)_\d+$/) { $refID=$1 }
    last if $refID eq $altID;
  }
print "TRACE ", __LINE__, "\n";
  $fastaWriter->writeFastaFromRef($refDef,$refSeqRef,$refFastaTemp);
  $fastaWriter->writeFastaFromRef($altDef,$altSeqRef,$altFastaTemp);
print "TRACE ", __LINE__, "\n";
  if($strand eq "-") {
    System("revcomp-fasta.pl $refFastaTemp > $tempRevcomp ; mv $tempRevcomp $refFastaTemp");
    System("revcomp-fasta.pl $altFastaTemp > $tempRevcomp ; mv $tempRevcomp $altFastaTemp");
  }
print "TRACE ", __LINE__, "\n";
  my $numTranscripts=@$transcripts;
  for(my $i=0 ; $i<$numTranscripts ; ++$i) {
    my $transcript=$transcripts->[$i];
    my $transID=$transcript->getTranscriptId();
    print "transcript $transID\n";
    next unless $STOP_AFTER eq "" || $STOP_AFTER eq $transID;
    ++$totalTranscripts;
    if($MAX_TRANSCRIPTS>0 && $totalTranscripts>$MAX_TRANSCRIPTS) { exit }
    my $gff=$transcript->toGff();
print "TRACE ", __LINE__, "\n";
    open(GFF,">$oneGeneGFF") || die "can't write file $oneGeneGFF";
    print GFF "$gff";
    close(GFF);
print "TRACE ", __LINE__, "\n";
    my $revFlag="";
    if($strand eq "-") {
      my $substrateLength=length($$refSeqRef);
print "TRACE ", __LINE__, "\n";
      System("revcomp-gff.pl $oneGeneGFF $substrateLength > $tempRevcomp ; mv $tempRevcomp $oneGeneGFF");
      $revFlag="-c";
print "TRACE ", __LINE__, "\n";
    }
    my $errorsFlag=$MAX_VCF_ERRORS>=0 ? "-e $MAX_VCF_ERRORS" : "";
    my $command="$ACE/ace $QUIET $errorsFlag $revFlag $modelFile $oneGeneGFF $refFastaTemp $altFastaTemp $gffTemp $aceTemp";
    print "$command\n";
    my $err=`$command`;
    #print "$err\n"
print "TRACE ", __LINE__, "\n";
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
      print "exiting for debugging\n";
      exit;
    }
print "TRACE ", __LINE__, "\n";
  }
print "TRACE ", __LINE__, "\n";
}
print "TRACE ", __LINE__, "\n";
unlink($refFastaTemp); unlink($altFastaTemp); 
#unlink($labelingTemp);
unlink($aceTemp); unlink($oneGeneGFF); unlink($tempRevcomp); unlink($gffTemp);
#unlink($tempXML);
#close(LABELING);
print "TRACE ", __LINE__, "\n";
close(ACE);
system("date");
print "ACE done.\n";

sub System
{
  my ($cmd)=@_;
  print("$cmd\n");
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
  #print "GC%=$rounded #GC=$GC #ACGT=$ATCG\n"; ###
  my $range;
  if($gc<=0.43) { $range="0-43" }
  elsif($gc<=0.51) { $range="43-51" }
  elsif($gc<=0.57) { $range="51-57" }
  else { $range="57-100" }
  return "$dir/ace.$range.config";
}






