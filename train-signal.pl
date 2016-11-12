#!/bin/env perl
use strict;
use TempFilename;
use FastaReader;
use FastaWriter;
use ProgramName;
use Getopt::Std;
our $opt_c;
getopts('c:');

my $STANDARD_CONTEXT_LENGTH=80; # what the example-generator script provides

# PROCESS COMMAND LINE
my $name=ProgramName::get();
my $usage=
"$name <model-type> <pos.fasta> <neg.fasta> <filestem>
     <\% train> <signal-type> <consensus-offset> 
     <consensus-length> <context-window-length> <min-sens> 
     <order> <min-sample-size> <boosting-iterations> 
     <boosting-percentile> [-c]
where
     <model-type> is one of {WMM,WAM,WWAM}
     <signal-type> is one of {ATG,TAG,GT,AG,PROMOTER,POLYA}
     <boosting-percentile> should be less than 0.2 for best results
     -c L = training file has L bases of context before and after 
            signals
";
die "$usage\n" unless @ARGV==14;
my ($modelType,$posFasta,$negFasta,$filestem,$percentTrain,$signalType,
    $consensusOffset,$consensusLength,$contextWindowLength,
    $minSensitivity,$order,$minSampleSize,$boostingIterations,
    $boostPercentile)=@ARGV;

# DO SOME INITIALIZATION
my $genezilla=$ENV{"GENEZILLA"};
if(!defined($genezilla)) 
  {die "environment variable GENEZILLA must be set to GeneZilla install directory\n"}
my $dashG=($percentTrain<1 ? "-g" : "");

my $posExamples=substringFasta($posFasta,$consensusOffset,
			       $contextWindowLength,"pos-$signalType");
my $negExamples=substringFasta($negFasta,$consensusOffset,
			       $contextWindowLength,"neg-$signalType");
my $dashB=($boostingIterations>0 ? "-b $boostingIterations" : "");
my $command="$genezilla/train-signal-sensor $dashG -o $order -s ".
  "$minSampleSize $modelType $posExamples ".
  "$negExamples $filestem $percentTrain $signalType ".
  "$consensusOffset $consensusLength $minSensitivity $dashB".
  " -p $boostPercentile";
print "$command\n";
system($command);
#unlink $posExamples;
#unlink $negExamples;


#---------------------------------------------------------------------
#---------------------------------------------------------------------
sub substringFasta
  {
    my ($infile,$consensusOffset,$contextWindowLength,$filestem)=@_;
    my $outfile="$filestem.tmp";#TempFilename::generate();
    print "reading \"$infile\"\n";
    my $reader=new FastaReader($infile);
    open(OUT,">$outfile") || die "can't create temp file $outfile";
    my $writer=new FastaWriter;
    my $begin=$STANDARD_CONTEXT_LENGTH-$consensusOffset;
    while(1)
      {
	my ($defline,$seq)=$reader->nextSequence();
	last unless defined $defline;
	my $subseq=substr($seq,$begin,$contextWindowLength);
	$writer->addToFasta($defline,$subseq,\*OUT);
      }
    close(OUT);
    return $outfile;
  }





