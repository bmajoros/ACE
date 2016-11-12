#!/bin/env perl
use strict;
use GffTranscriptReader;
use FastaReader;
use FastaWriter;

my $iso="0-100";

my $MAX_EXAMPLES=1000;#1000000;#2000;
my $DONOR_CONTEXT_BEGIN=-80;
my $DONOR_CONTEXT_LENGTH=162;
my $ACCEPTOR_CONTEXT_BEGIN=-80;
my $ACCEPTOR_CONTEXT_LENGTH=162;
my $START_CONTEXT_BEGIN=-80;
my $START_CONTEXT_LENGTH=163;
my $STOP_CONTEXT_BEGIN=-80;
my $STOP_CONTEXT_LENGTH=163;

my $numStartCodons=`grep -c '>' start-codons.fasta`;
my $numStopCodons=`grep -c '>' stop-codons.fasta`;
my $numDonors=`grep -c '>' donors.fasta`;
my $numAcceptors=`grep -c '>' acceptors.fasta`;

open(DONORS,">non-donors.fasta");
open(ACCEPTORS,">non-acceptors.fasta");
open(STARTCODONS,">non-start-codons.fasta");
open(STOPCODONS,">non-stop-codons.fasta");

my $exonHash=FastaReader::readAll("internal-exons.fasta");
my @exons=values %$exonHash;
my $intergenicHash=FastaReader::readAll("intergenic.fasta");
my @intergenic=values %$intergenicHash;

#######################################################################
# Non-Donors
#######################################################################
my $donors=
  getExampleSignals(\@exons,"GT",$DONOR_CONTEXT_BEGIN,$DONOR_CONTEXT_LENGTH,
	      $MAX_EXAMPLES/2);
my $moreDonors=
  getExampleSignals(\@intergenic,"GT",$DONOR_CONTEXT_BEGIN,
		    $DONOR_CONTEXT_LENGTH,$MAX_EXAMPLES/2);
push @$donors,@$moreDonors;
while(@$donors>$numDonors) {pop @$donors}
foreach my $donor (@$donors)
  {
    print DONORS ">non-donor\n$donor\n";
  }

#######################################################################
# Non-Acceptors
#######################################################################
my $acceptors=
  getExampleSignals(\@exons,"AG",$ACCEPTOR_CONTEXT_BEGIN,
		    $ACCEPTOR_CONTEXT_LENGTH,$MAX_EXAMPLES/2);
my $moreAcceptors=
  getExampleSignals(\@intergenic,"AG",$ACCEPTOR_CONTEXT_BEGIN,
	      $ACCEPTOR_CONTEXT_LENGTH,$MAX_EXAMPLES/2);
push @$acceptors,@$moreAcceptors;
while(@$acceptors>$numAcceptors) {pop @$acceptors}
foreach my $acceptor (@$acceptors)
  {
    print ACCEPTORS ">non-acceptor\n$acceptor\n";
  }

#######################################################################
# Non-Start codons
#######################################################################
my $starts=
  getExampleSignals(\@exons,"ATG",$START_CONTEXT_BEGIN,$START_CONTEXT_LENGTH,
	      $MAX_EXAMPLES/2);
my $moreStarts=
  getExampleSignals(\@intergenic,"ATG",$START_CONTEXT_BEGIN,
		    $START_CONTEXT_LENGTH,$MAX_EXAMPLES/2);
push @$starts,@$moreStarts;
while(@$starts>$numStartCodons) {pop @$starts}
foreach my $start (@$starts)
  {
    print STARTCODONS ">non-start-codon\n$start\n";
  }

#######################################################################
# Non-Stop codons
#######################################################################
my $stops=
  getExampleSignals(\@exons,"TAG",$STOP_CONTEXT_BEGIN,$STOP_CONTEXT_LENGTH,
	      $MAX_EXAMPLES/6);
my $moreStops=
  getExampleSignals(\@exons,"TGA",$STOP_CONTEXT_BEGIN,$STOP_CONTEXT_LENGTH,
	      $MAX_EXAMPLES/6);
push @$stops,@$moreStops;
$moreStops=
  getExampleSignals(\@exons,"TAA",$STOP_CONTEXT_BEGIN,$STOP_CONTEXT_LENGTH,
	      $MAX_EXAMPLES/6);
push @$stops,@$moreStops;
$moreStops=
  getExampleSignals(\@intergenic,"TAG",$STOP_CONTEXT_BEGIN,
		    $STOP_CONTEXT_LENGTH,$MAX_EXAMPLES/6);
push @$stops,@$moreStops;
$moreStops=
  getExampleSignals(\@intergenic,"TGA",$STOP_CONTEXT_BEGIN,
		    $STOP_CONTEXT_LENGTH,$MAX_EXAMPLES/6);
push @$stops,@$moreStops;
$moreStops=
  getExampleSignals(\@intergenic,"TAA",$STOP_CONTEXT_BEGIN,
		    $STOP_CONTEXT_LENGTH,$MAX_EXAMPLES/6);
push @$stops,@$moreStops;
while(@$stops>$numStopCodons) {pop @$stops}
foreach my $stop (@$stops)
  {
    print STOPCODONS ">non-stop-codon\n$stop\n";
  }

close(DONORS);
close(ACCEPTORS);
close(STARTCODONS);
close(STOPCODONS);

#######################################################################
# Non-Exons
#######################################################################
my $pool;
appendToPool("introns.fasta",\$pool);
appendToPool("intergenic.fasta",\$pool);
getExampleContent("initial-exons",\$pool);
getExampleContent("internal-exons",\$pool);
getExampleContent("final-exons",\$pool);
getExampleContent("single-exons",\$pool);

#######################################################################
# Non-Intergenic
#######################################################################
undef $pool;
appendToPool("initial-exons.fasta",\$pool);
appendToPool("final-exons.fasta",\$pool);
appendToPool("single-exons.fasta",\$pool);
getExampleContent("intergenic",\$pool);

#######################################################################
# Non-Introns
#######################################################################
undef $pool;
appendToPool("internal-exons.fasta",\$pool);
appendToPool("initial-exons.fasta",\$pool);
appendToPool("final-exons.fasta",\$pool);
getExampleContent("introns",\$pool);

#------------------------------------------------------------
sub appendToPool
  {
    my ($infile,$pool)=@_;
    my $reader=new FastaReader($infile);
    while(1)
      {
	my ($defline,$sequence)=$reader->nextSequence();
	if(!defined($defline)) {last}
	if($sequence=~/\n/) {chop $sequence}
	$$pool.=$sequence;
      }
  }
#------------------------------------------------------------
sub getExampleContent
  {
    my ($type,$pool)=@_;
    my $nonType="non-$type";
    my $poolLen=length($$pool);
    my $writer=new FastaWriter;
    open(OUT,">$nonType.fasta") || die;
    my $pos=0;
    my $i=0;
    my $reader=new FastaReader("$type.fasta");
    while(1)
      {
	my ($defline,$sequence)=$reader->nextSequence();
	if(!defined($defline)) {last}
	my $len=length($sequence);
	if($pos+$len>=$poolLen) {$pos=int(rand(100))}
	my $seq=substr($$pool,$pos,$len);
	$pos+=$len;
	$writer->addToFasta(">$nonType (G+C: $iso) #$i",$seq,\*OUT);
	++$i;
      }
    close(OUT);
  }
#------------------------------------------------------------
sub getExampleSignals
  {
    my ($source,$signal,$contextBegin,$contextLength,$max)=@_;
    my $signalLen=length($signal);
    my $left=-$contextBegin;
    my $right=$contextLength-$left-$signalLen;
    my @examples;
    my $n=@$source;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $exon=$source->[$i];
	while(1)
	  {
	    if($exon=~/(.{$left}$signal.{$right})(.*)/)
	      {
		my ($example,$remainder)=($1,$2);
		push @examples,$example;
		$exon=$remainder;
		next;
	      }
	    last;
	  }
	last unless @examples<$max;
      }
    return \@examples;
  }
#--------------------------------------------------------------
