#!/bin/env perl
use strict;
use GffTranscriptReader;
use FastaReader;
use FastaWriter;
use ProgramName;
$|=1;

my $MAX_DUPLICATE_COUNT=100; # set to 1 to eliminate duplicates

# PROCESS THE COMMAND LINE
my $name=ProgramName::get();
my $usage="$name <*.gff> <*.fasta> <GT,GC,AT> <AG,AC> <ATG> <TGA,TAA,TAG>";
die "$usage\n" unless @ARGV==6;
my ($gffFilename,$fastaFilename,$donors,$acceptors,$startCodons,$stopCodons)
  =@ARGV;

my @donors=split/,/,$donors;
my @acceptors=split/,/,$acceptors;
my @startCodons=split/,/,$startCodons;
my @stopCodons=split/,/,$stopCodons;
my (%donors,%acceptors,%startCodons,%stopCodons);
foreach my $x (@donors) {$donors{$x}=1}
foreach my $x (@acceptors) {$acceptors{$x}=1}
foreach my $x (@startCodons) {$startCodons{$x}=1}
foreach my $x (@stopCodons) {$stopCodons{$x}=1}

my $shouldTrim=0;
#my @stopCodons=("TGA","TAA","TAG");
#my %stopCodons;
#foreach my $codon (@stopCodons) {$stopCodons{$codon}=1}

# SET SOME LIMITS
my $MAX_EXONS=1000000; #10000;
my $MAX_DONORS=1000000; #500000;
my $MAX_ACCEPTORS=1000000; #500000;
my $MAX_START_CODONS=100000; #50000;
my $MAX_STOP_CODONS=1000; #50000;
my $MAX_INTRONS=1000000; #10000;
my $MAX_INTERGENICS=1000; #10000;
my $MAX_UTRS=1000; #10000;
my $MAX_TOTAL_LENGTH=30000000; #1000000;
my $DONOR_CONTEXT_BEGIN=-80;
my $DONOR_CONTEXT_LENGTH=162;
my $ACCEPTOR_CONTEXT_BEGIN=-80;
my $ACCEPTOR_CONTEXT_LENGTH=162;
my $START_CONTEXT_BEGIN=-80;
my $START_CONTEXT_LENGTH=163;
my $STOP_CONTEXT_BEGIN=-80;
my $STOP_CONTEXT_LENGTH=163;
my $MAX_INTERGENIC=10000;
my $UTR_EXTENT=50;
my ($numInitialExons,$numFinalExons,$numInternalExons,$numSingleExons,
    $numDonors,$numAcceptors,$numStartCodons,$numStopCodons,
    $numIntrons,$numIntergenic,$totalInitialExon,$totalInternalExon,
    $totalFinalExon,$totalSingleExon,$totalIntron,$totalIntergenic,
    $numFivePrimeUTR,$numThreePrimeUTR);

# LOAD ALL THE TRANSCRIPTS IN THE GFF FILE
my ($isoLowerBound,$isoUpperBound)=(0,100);
my $iso="$isoLowerBound\-$isoUpperBound";
my $gffReader=new GffTranscriptReader();
$gffReader->setStopCodons(\%stopCodons);
my $transcripts=$gffReader->loadGFF($gffFilename);
my $writer=new FastaWriter;

# CREATE OUTPUT FILES
open(INITIALEXONS,">initial-exons.fasta");
open(FINALEXONS,">final-exons.fasta");
open(INTERNALEXONS,">internal-exons.fasta");
open(SINGLEEXONS,">single-exons.fasta");
open(DONORS,">donors.fasta");
open(ACCEPTORS,">acceptors.fasta");
open(STARTCODONS,">start-codons.fasta");
open(STOPCODONS,">stop-codons.fasta");
open(INTRONS,">introns.fasta");
open(INTERGENIC,">intergenic.fasta");
open(FIVEPRIMEUTR,">five-prime-UTR.fasta");
open(THREEPRIMEUTR,">three-prime-UTR.fasta");

# PUT ALL THE TRANSCRIPTS IN A HASH TABLE INDEXED BY CONTIG ID
my %transcriptsOnContig;
my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $contig=$transcript->getSubstrate();
    push @{$transcriptsOnContig{$contig}},$transcript;
  }

# NOW READ THE CONTIGS ONE BY ONE, PROCESSING ALL THE TRANSCRIPTS ON
# EACH CONTIG AS WE READ IT
my $fastaReader=new FastaReader($fastaFilename);
while(1)
  {
    my ($defline,$substrateSeq)=$fastaReader->nextSequence();
    last unless defined $defline;
    $defline=~/>(\S+)/ || die "Can't parse defline: $defline\n";
    my $substrateId=$1;
    my $transcripts=$transcriptsOnContig{$substrateId};
    next unless defined $transcripts;
    my $numTranscripts=@$transcripts;
    my $substrateLen=length $substrateSeq;
    my %featuresSeen; # coordinate strings of the form: "$begin,$end"
    my %featureTypes;
    for(my $i=0 ; $i<$numTranscripts ; ++$i)
      {
	my $transcript=$transcripts->[$i];
	next unless transcriptIsOk($transcript);
	my $substrate=$transcript->{substrate};
	if(!defined $substrateSeq) 
	  {print "SKIPPING $substrate DUE TO SUBSTRATE\n";next;}
	my $transcriptId=$transcript->{transcriptId};
	my $strand=$transcript->{strand};
	
	if($shouldTrim)
	  {
	    $transcript->trimUTR(\$substrateSeq);
	    $transcript->recomputeBoundaries();
	  }
	
	##################################################################
	# Exons
	##################################################################
	if(!$transcript->loadExonSequences(\$substrateSeq))
	  {print "SKIPPING DUE TO SEQUENCE PROBLEMS\n";next;}
	my $numExons=$transcript->numExons();
	if($numExons==1)
	  {
	    my $exon=$transcript->getIthExon(0);
	    my $sequence=$exon->getSequence();
	    my $substrLen=$exon->getLength()-6;
	    my $sigBegin=$exon->getBegin()+3;
	    my $signature="$sigBegin,$substrLen";
	    $sequence=substr($sequence,3,$substrLen);
	    my $length=length($sequence);
	    $featureTypes{$signature}="single-exon";
	    if($numSingleExons<$MAX_EXONS && 
	       $totalSingleExon<$MAX_TOTAL_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $sequence!~/N/ &&
	       $exon->getType() eq "single-exon")
	      {
		$writer->addToFasta(
			       ">$transcriptId /frame=0 /length=$length",
				    $sequence,\*SINGLEEXONS);
		++$numSingleExons;
		$totalSingleExon+=$length;
	      }
	  }
	else
	  {
	    my $firstExon=$transcript->getIthExon(0);
	    my $sequence=$firstExon->getSequence();
	    my $substrLen=$firstExon->getLength()-3;
	    my $sigBegin=$firstExon->getBegin()+3;
	    my $signature="$sigBegin,$substrLen";
	    $sequence=substr($sequence,3,$substrLen);
	    my $length=length($sequence);
	    $featureTypes{$signature}="initial-exon";
	    if($numInitialExons<$MAX_EXONS &&
	       $totalInitialExon<$MAX_TOTAL_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $sequence!~/N/ &&
	       $firstExon->getType() eq "initial-exon")
	      {
		$writer->addToFasta(
		     ">$transcriptId /exon=0 /frame=0 /length=$length",
				    $sequence,\*INITIALEXONS);
		++$numInitialExons;
		$totalInitialExon+=$length;
	      }
	    
	    my $exonNum=$numExons-1;
	    my $lastExon=$transcript->getIthExon($exonNum);
	    my $frame=$lastExon->{frame};
	    $sequence=$lastExon->getSequence();
	    my $substrLen=$lastExon->getLength()-3;
	    my $sigBegin=$lastExon->getBegin();
	    my $signature="$sigBegin,$substrLen";
	    $sequence=substr($sequence,0,$substrLen);
	    my $length=length($sequence);
	    $featureTypes{$signature}="final-exon";
	    if($numFinalExons<$MAX_EXONS &&
	       $totalFinalExon<$MAX_TOTAL_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $sequence!~/N/ &&
	       $lastExon->getType() eq "final-exon")
	      {
		$writer->addToFasta(
	     ">$transcriptId /exon=$exonNum /frame=$frame /length=$length",
				    $sequence,\*FINALEXONS);
		++$numFinalExons;
		$totalFinalExon+=$length;
	      }
	    
	    my $numInternal=$numExons-2;
	    for(my $j=0 ; $j<$numInternal ; ++$j)
	      {
		$exonNum=$j+1;
		my $exon=$transcript->getIthExon($exonNum);
		my $frame=$exon->{frame};
		$sequence=$exon->getSequence();
		my $length=length($sequence);
		my $sigBegin=$lastExon->getBegin();
		my $signature="$sigBegin,$length"; #"0,$length";
		$featureTypes{$signature}="internal-exon";
		if($numInternalExons<$MAX_EXONS &&
		   $totalInternalExon<$MAX_TOTAL_LENGTH &&
		   ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $sequence!~/N/)
		  {
		    $writer->addToFasta(
              ">$transcriptId /exon=$exonNum /frame=$frame /length=$length",
					$sequence,\*INTERNALEXONS);
		    ++$numInternalExons;
		    $totalInternalExon+=$length;
		  }
	      }
	  }
	
	####################################################################
	# Splice sites
	####################################################################
	if($numExons>1)
	  {
	    my $firstExon=$transcript->getIthExon(0);
	    my $sigBegin=$firstExon->getEnd();
	    my $signature="$sigBegin,2";
	    my $context=getDonor($firstExon,\$substrateSeq,$strand);
	    $featureTypes{$signature}="donor";

	    #####
	    my $contextLen=length($context);
	    my $debugFeaturesSeen=$featuresSeen{$signature};
	    my $containsNs=($context=~/N/);
	    #####

	    my $consensus=substr($context,80,2);
#	    if(($consensus eq "GT" || $consensus eq "GC") &&
	    if($donors{$consensus} &&
	       $numDonors<$MAX_DONORS &&
	       length($context)==$DONOR_CONTEXT_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $context!~/N/)
	      {
		$writer->addToFasta(">$transcriptId",$context,\*DONORS);
		++$numDonors;
	      }

	    my $exonNum=$numExons-1;
	    my $lastExon=$transcript->getIthExon($exonNum);
	    $sigBegin=$lastExon->getBegin()-2;
	    my $signature="$sigBegin,2";
	    $context=getAcceptor($lastExon,\$substrateSeq,$strand);
	    $featureTypes{$signature}="acceptor";

	    my $consensus=substr($context,80,2);
#	    if($consensus eq "AG" && 
	    if($acceptors{$consensus} &&
	       $numAcceptors<$MAX_ACCEPTORS &&
	       length($context)==$ACCEPTOR_CONTEXT_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $context!~/N/)
	      {
		$writer->addToFasta(">$transcriptId",$context,\*ACCEPTORS);
		++$numAcceptors;
	      }
	    
	    my $numInternal=$numExons-2;
	    for(my $j=0 ; $j<$numInternal ; ++$j)
	      {
		$exonNum=$j+1;
		my $exon=$transcript->getIthExon($exonNum);
		$sigBegin=$exon->getEnd();
		$signature="$sigBegin,2";
		$context=getDonor($exon,\$substrateSeq,$strand);
		$featureTypes{$signature}="donor";

	    #####
	    my $contextLen=length($context);
	    my $debugFeaturesSeen=$featuresSeen{$signature};
	    my $containsNs=($context=~/N/);
	    #####

		my $consensus=substr($context,80,2);
#		if(($consensus eq "GT" || $consensus eq "GC") &&
		if($donors{$consensus} &&
		   $numDonors<$MAX_DONORS &&
		   length($context)==$DONOR_CONTEXT_LENGTH &&
		   ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
		   $context!~/N/)
		  {
		    $writer->addToFasta(">$transcriptId",$context,\*DONORS);
		    ++$numDonors;
		  }
		$context=getAcceptor($exon,\$substrateSeq,$strand);
		$sigBegin=$exon->getBegin()-2;
		$signature="$sigBegin,2";
		$featureTypes{$signature}="acceptor";

		my $consensus=substr($context,80,2);
#		if($consensus eq "AG" &&
		if($acceptors{$consensus} && 
		   $numAcceptors<$MAX_ACCEPTORS &&
		   length($context)==$ACCEPTOR_CONTEXT_LENGTH &&
		   ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
		   $context!~/N/)
		  {
		    $writer->addToFasta(">$transcriptId",$context,
					\*ACCEPTORS);
		    ++$numAcceptors;
		  }
	      }
	  }
	
	###################################################################
	# Start & stop codons
	###################################################################
	my $firstExon=$transcript->getIthExon(0);
	my $lastExon=$transcript->getIthExon($numExons-1);
	if($strand eq "+")
	  {
	    my $featureCoord=$firstExon->{begin};
	    my $begin=$featureCoord+$START_CONTEXT_BEGIN;
	    my $signature="$begin,3";
	    my $context=substr($substrateSeq,$begin,$START_CONTEXT_LENGTH);
	    my $consensus=substr($substrateSeq,$featureCoord,3);
	    $featureTypes{$signature}="start-codon";
	    if($numStartCodons<$MAX_START_CODONS && 
	       length($context)==$START_CONTEXT_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $startCodons{$consensus} &&
	       $context!~/N/ &&
	       ($firstExon->getType() eq "initial-exon" ||
	       $firstExon->getType() eq "single-exon"))
	      {
		$writer->addToFasta(">$transcriptId",$context,
				    \*STARTCODONS);
		++$numStartCodons;
	      }

	    $featureCoord=$lastExon->{end}-3;
	    $begin=$featureCoord+$STOP_CONTEXT_BEGIN;
	    $signature="$begin,3";
	    $context=substr($substrateSeq,$begin,$STOP_CONTEXT_LENGTH);
	    my $consensus=substr($substrateSeq,$featureCoord,3);
	    $featureTypes{$signature}="stop-codon";
	    if($numStopCodons<$MAX_STOP_CODONS &&
	       length($context)==$STOP_CONTEXT_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $stopCodons{$consensus} &&
	       $context!~/N/ &&
	       ($lastExon->getType() eq "final-exon" ||
		$lastExon->getType() eq "single-exon"))
	      {
		$writer->addToFasta(">$transcriptId",$context,\*STOPCODONS);
		++$numStopCodons;
	      }
	  }
	else
	  {
	    my $featureEnd=$firstExon->{end};
	    my $begin=
	      $featureEnd-$START_CONTEXT_BEGIN-$START_CONTEXT_LENGTH;
	    my $signature="$begin,3";
	    my $context=substr($substrateSeq,$begin,$START_CONTEXT_LENGTH);
	    $context=Translation::reverseComplement(\$context);
	    $featureTypes{$signature}="start-codon";
	    my $consensus=substr($substrateSeq,$begin-$START_CONTEXT_BEGIN,3);
	    $consensus=Translation::reverseComplement(\$consensus);
	    if($numStartCodons<$MAX_START_CODONS && 
	       length($context)==$START_CONTEXT_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $startCodons{$consensus} &&
	       $context!~/N/ &&
	       ($firstExon->getType() eq "initial-exon" ||
		$firstExon->getType() eq "single-exon"))
	      {
		$writer->addToFasta(">$transcriptId",$context,
				    \*STARTCODONS);
		++$numStartCodons;
	      }
	    
	    $featureEnd=$lastExon->{begin}+3;
	    $begin=$featureEnd-$STOP_CONTEXT_BEGIN-$STOP_CONTEXT_LENGTH;
	    $signature="$begin,3";
	    $context=substr($substrateSeq,$begin,$STOP_CONTEXT_LENGTH);
	    $context=Translation::reverseComplement(\$context);
	    my $consensus=substr($substrateSeq,$begin-$STOP_CONTEXT_BEGIN,3);
	    $consensus=Translation::reverseComplement(\$consensus);
	    $featureTypes{$signature}="stop-codon";
	    if($numStopCodons<$MAX_STOP_CODONS &&
	       length($context)==$STOP_CONTEXT_LENGTH &&
	       ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
	       $stopCodons{$consensus} && 
	       $context!~/N/ &&
	       ($lastExon->getType() eq "final-exon" ||
		$lastExon->getType() eq "single-exon"))
	      {
		$writer->addToFasta(">$transcriptId",$context,\*STOPCODONS);
		++$numStopCodons;
	      }
	  }
	
	####################################################################
	# Introns
	####################################################################
	if($numExons>1)
	  {
	    for(my $j=0 ; $j<$numExons-1 ; ++$j)
	      {
		my $thisExon=$transcript->getIthExon($j);
		my $nextExon=$transcript->getIthExon($j+1);
		my ($intron,$sigBegin,$sigLen);
		if($strand eq "+")
		  {
		    $sigBegin=$thisExon->{end};
		    $sigLen=$nextExon->{begin}-$thisExon->{end};
		    $intron=substr($substrateSeq,$sigBegin,$sigLen);
		  }
		else
		  {
		    $sigBegin=$nextExon->{end};
		    $sigLen=$thisExon->{begin}-$nextExon->{end};
		    $intron=substr($substrateSeq,$sigBegin,$sigLen);
		    $intron=Translation::reverseComplement(\$intron);
		  }
		my $signature="$sigBegin,$sigLen";
		$featureTypes{$signature}="intron";
		if($numIntrons<$MAX_INTRONS &&
		   $totalIntron<$MAX_TOTAL_LENGTH &&
		   ++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
		   $intron!~/N/)
		  {
		    my $length=length($intron);
		    $writer->addToFasta(
                          ">$transcriptId /intron=$j /length=$length"
					,$intron,
					\*INTRONS);
		    ++$numIntrons;
		    $totalIntron+=$length;
		  }
	      }
	  }
	
	####################################################################
	# UTR
	####################################################################
	if($strand eq "+")
	  {
	    if($numFivePrimeUTR<$MAX_UTRS)
	      {
		my $utrBegin=$firstExon->{begin}-$UTR_EXTENT;
		if($utrBegin<0) {$utrBegin=0};
		my $utrExtent=$firstExon->{begin}-$utrBegin;
		my $signature="$utrBegin,$utrExtent";
		$featureTypes{$signature}="5'-UTR";
		my $fivePrime=
		  substr($substrateSeq,$utrBegin,$utrExtent);
		if(++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
		   $fivePrime!~/N/ &&
		   ($firstExon->getType() eq "initial-exon" ||
		    $firstExon->getType() eq "single-exon"))
		  {
		    $writer->addToFasta(">$transcriptId",
					$fivePrime,\*FIVEPRIMEUTR);
		    ++$numFivePrimeUTR;
		  }
	      }
	    if($numThreePrimeUTR<$MAX_UTRS)
	      {
		my $utrBegin=$lastExon->{end}+1;
		my $utrExtent=$substrateLen-$utrBegin;
		if($utrExtent>$UTR_EXTENT) {$utrExtent=$UTR_EXTENT}
		my $signature="$utrBegin,$utrExtent";
		$featureTypes{$signature}="3'-UTR";
		my $threePrime=
		  substr($substrateSeq,$utrBegin,$utrExtent);
		if(++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
		   $threePrime!~/N/ &&
		   ($lastExon->getType() eq "final-exon" ||
		    $lastExon->getType() eq "single-exon"))
		  {
		    $writer->addToFasta(">$transcriptId",$threePrime,
					\*THREEPRIMEUTR);
		    ++$numThreePrimeUTR;
		  }
	      }
	  }
	else # REVERSE STRAND:
	  {
	    if($numFivePrimeUTR<$MAX_UTRS)
	      {
		my $utrBegin=$firstExon->{end}+1;
		my $utrExtent=$substrateLen-$utrBegin;
		if($utrExtent>$UTR_EXTENT) {$utrExtent=$UTR_EXTENT}
		my $signature="$utrBegin,$utrExtent";
		$featureTypes{$signature}="5'-UTR";
		my $fivePrime=
		  substr($substrateSeq,$utrBegin,$utrExtent);
		if(++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
		   $fivePrime!~/N/ &&
		   ($firstExon->getType() eq "initial-exon" ||
		    $firstExon->getType() eq "single-exon"))
		  {
		    $fivePrime=Translation::reverseComplement(\$fivePrime);
		    $writer->addToFasta(">$transcriptId",$fivePrime,
					\*FIVEPRIMEUTR);
		    ++$numFivePrimeUTR;
		  }
	      }
	    if($numThreePrimeUTR<$MAX_UTRS)
	      {
		my $utrBegin=$lastExon->{begin}-$UTR_EXTENT;
		if($utrBegin<0) {$utrBegin=0};
		my $utrExtent=$lastExon->{begin}-$utrBegin;
		my $signature="$utrBegin,$utrExtent";
		$featureTypes{$signature}="3'-UTR";
		my $threePrime=
		  substr($substrateSeq,$utrBegin,$utrExtent);
		if(++$featuresSeen{$signature}<=$MAX_DUPLICATE_COUNT &&
		   $threePrime!~/N/ &&
		   ($lastExon->getType() eq "final-exon" ||
		    $lastExon->getType() eq "single-exon"))
		  {
		    $threePrime=
		      Translation::reverseComplement(\$threePrime);
		    $writer->addToFasta(">$transcriptId",$threePrime,
					\*THREEPRIMEUTR);
		    ++$numThreePrimeUTR;
		  }
	      }
	  }
	
	####################################################################
	# INTERGENIC
	####################################################################
	if($numIntergenic<$MAX_INTERGENICS &&
	   $totalIntergenic<$MAX_TOTAL_LENGTH)
	  {
	    my $endIntergenic=$transcript->getBegin()-$UTR_EXTENT;
	    my $beginIntergenic=$endIntergenic-$MAX_INTERGENIC;
	    if($beginIntergenic<0) {$beginIntergenic=0}
	    my $intergenicLen=$endIntergenic-$beginIntergenic;
	    my $intergenicSeq=
	      substr($substrateSeq,$beginIntergenic,$intergenicLen);
	    if($intergenicLen>0 &&
	       $intergenicSeq!~/N/ &&
	       ($transcript->getIthExon(0)->getType() eq "initial-exon" ||
		$transcript->getIthExon(0)->getType() eq "single-exon"))
	      {
		$writer->addToFasta(">$transcriptId",
				    $intergenicSeq,\*INTERGENIC);
		++$numIntergenic;
		$totalIntergenic+=$intergenicLen;
	      }

	    $beginIntergenic=$transcript->getEnd()+$UTR_EXTENT;
	    $endIntergenic=$beginIntergenic+$MAX_INTERGENIC;
	    if($endIntergenic>$substrateLen) {$endIntergenic=$substrateLen}
	    $intergenicLen=$endIntergenic-$beginIntergenic;
	    my $intergenicSeq=
	      substr($substrateSeq,$beginIntergenic,$intergenicLen);
	    if($intergenicLen>0 &&
	       $intergenicSeq!~/N/ &&
	       ($transcript->getIthExon($transcript->numExons()-1)->
		getType() eq "final-exon" ||
		$transcript->getIthExon($transcript->numExons()-1)->
		getType() eq "single-exon"))
	      {
		$writer->addToFasta(">$transcriptId",
				    $intergenicSeq,\*INTERGENIC);
		++$numIntergenic;
		$totalIntergenic+=$intergenicLen;
	      }
	  }
      }

    my @signatures=keys %featuresSeen;
    foreach my $signature (@signatures)
      {
	my $count=$featuresSeen{$signature};
	if($count>1)
	  {
	    my $featureType=$featureTypes{$signature};
	    #print "$featureType ($signature) $count\n";
	  }
      }
  }

close(INITIALEXONS);
close(FINALEXONS);
close(INTERNALEXONS);
close(SINGLEEXONS);
close(DONORS);
close(ACCEPTORS);
close(STARTCODONS);
close(STOPCODONS);
close(INTRONS);
close(INTERGENIC);

#--------------------------------------------------------------------------
sub getDonor
  {
    my ($exon,$chunk,$strand)=@_;
    if($strand eq "+")
      {
	my $featureCoord=$exon->{end};
	my $begin=$featureCoord+$DONOR_CONTEXT_BEGIN;
	my $context=substr($$chunk,$begin,$DONOR_CONTEXT_LENGTH);
	return $context;
      }
    else
      {
	my $featureEnd=$exon->{begin};
	my $begin=$featureEnd-$DONOR_CONTEXT_BEGIN-$DONOR_CONTEXT_LENGTH;
	my $context=substr($$chunk,$begin,$DONOR_CONTEXT_LENGTH);
	$context=Translation::reverseComplement(\$context);
	return $context;
      }
  }
#--------------------------------------------------------------------------
sub getAcceptor
  {
    my ($exon,$chunk,$strand)=@_;
    if($strand eq "+")
      {
	my $featureCoord=$exon->{begin}-2;
	my $begin=$featureCoord+$ACCEPTOR_CONTEXT_BEGIN;
	my $context=substr($$chunk,$begin,$ACCEPTOR_CONTEXT_LENGTH);
	return $context;
      }
    else
      {
	my $featureEnd=$exon->{end}+2;
	my $begin=$featureEnd-$ACCEPTOR_CONTEXT_BEGIN-
	  $ACCEPTOR_CONTEXT_LENGTH;
	my $context=substr($$chunk,$begin,$ACCEPTOR_CONTEXT_LENGTH);
	$context=Translation::reverseComplement(\$context);
	return $context;
      }
  }
#--------------------------------------------------------------------------
sub transcriptIsOk
  {
    my ($transcript)=@_;
    my $n=$transcript->numExons();
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $exon=$transcript->getIthExon($i);
	my $len=$exon->getEnd()-$exon->getBegin();
	if($len<3) {return 0}
      }
    return 1;
  }
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
