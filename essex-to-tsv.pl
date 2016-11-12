#!/bin/env perl
$|=1;
use strict;
use EssexParser;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.essex>\n" unless @ARGV==1;
my ($infile)=@ARGV;

print "GeneID\tTranscriptID\tRefSubstrate\tAltSubstrate\tVcfWarnings\tVcfErrors\tBegin\tEnd\tStrand\tType\tMapped\tSplicingChanges\tBrokenSpliceSite\tProteinDiffers\tFrameshift\tPrematureStop\tNMD\tStartCodonChange\tNoStartCodon\tNoTranscript\tNonstopDecay\n";

my $parser=new EssexParser($infile);
while(1) {
  my $report=$parser->nextElem();
  last unless $report;
  my $status=$report->findChild("status");
  next unless $status;
  my $code=$status->getIthElem(0);
  next unless $code;
  my $refTrans=$report->findChild("reference-transcript");
  next unless $refTrans;

  my $substrate=$report->getAttribute("substrate");
  my $transcriptID=$report->getAttribute("transcript-ID");
  my $geneID=$report->getAttribute("gene-ID");
  my $vcfWarnings=0+$report->getAttribute("vcf-warnings");
  my $vcfErrors=0+$report->getAttribute("vcf-effors");

  my $type=$refTrans->getAttribute("type");
  my $refSubstrate=$refTrans->getAttribute("substrate");
  my $strand=$refTrans->getAttribute("strand");
  my $begin=$refTrans->getAttribute("begin");
  my $end=$refTrans->getAttribute("end");

  my $splicingChanges=$status->hasDescendentOrDatum("splicing-changes");
  my $brokenSpliceSite=$status->hasDescendentOrDatum("broken-donor") ||
    $status->hasDescendentOrDatum("broken-acceptor");
  my $proteinDiffers=$status->hasDescendentOrDatum("protein-differs");
  my $frameshift=$status->hasDescendentOrDatum("frameshift");
  my $prematureStop=$status->hasDescendentOrDatum("premature-stop");
  my $NMD=$status->hasDescendentOrDatum("NMD");
  my $startCodonChange=$status->hasDescendentOrDatum("start-codon-change");
  my $noTranscript=$status->hasDescendentOrDatum("no-transcript");
  my $noStartCodon=$status->hasDescendentOrDatum("no-start-codon");
  my $mapped=$status->hasDescendentOrDatum("mapped");
  my $nonstopDecay=$status->hasDescendentOrDatum("nonstop-decay");

  my $mappedTrans=$report->findChild("mapped-transcript");
  my $lossOfCoding=$mappedTrans && $mappedTrans->getAttribute("type") ne $type
    ? 1 : 0;

  print "$geneID\t$transcriptID\t$refSubstrate\t$substrate\t$vcfWarnings\t$vcfErrors\t$begin\t$end\t$strand\t$type\t$mapped\t$splicingChanges\t$brokenSpliceSite\t$proteinDiffers\t$frameshift\t$prematureStop\t$NMD\t$startCodonChange\t$noStartCodon\t$noTranscript\t$nonstopDecay\n";

}



