#!/bin/env perl
use strict;
use EssexParser;
use ProgramName;

my $VERBOSE=0;

my $name=ProgramName::get();
die "$name <in.essex>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my (%counts,%flags);
my %NMD; my %noTranscriptBrokenSS; my %noTranscriptWeakenedSS;
my $brokenDonor=0; my $brokenAcceptor=0; my $weakenedDonor=0;
my $crypticSites=0; my $crypticSitesSquared=0; my $singleBrokenSites=0;
my $exonSkipping=0; my $exonSkippingSquared=0;
my $noStartCodon=0;
my $weakenedAcceptor=0; my %splicingChangesNMD; my %splicingChangesCoding;
my %splicingChangesNoncoding; my @altTranscripts; my @altSquared;
my %proteinDiffers; my %frameshift; my %nonstopDecay;
my %prematureStop; my %startCodonChange;
my @frameshifts; my @squaredFrameshifts;
my $parser=new EssexParser($infile);
while(1) {
  my $report=$parser->nextElem();
  last unless $report;
  my $status=$report->findChild("status");
  next unless $status;
  my $code=$status->getIthElem(0);
  next unless $code;
  my $gene=$report->getAttribute("gene-ID");
  my $coding;
  my $ref=$report->findChild("reference-transcript") || die "no reference";
  my $type=$ref->getAttribute("type");
  if($type eq "") { die "no type" }
  if($type="protein-coding") { $coding=1 }
  $brokenDonor+=$status->countDescendents("broken-donor");
  $brokenAcceptor+=$status->countDescendents("broken-acceptor");
  $weakenedDonor+=$status->countDescendents("weakened-donor");
  $weakenedAcceptor+=$status->countDescendents("weakened-acceptor");
  if($status->hasDescendentOrDatum("no-start-codon")) { ++$noStartCodon }
  my $altTranscripts=$status->findDescendents("transcript");
  if($altTranscripts && @$altTranscripts) {
    my $n=0+@$altTranscripts;
    push @altTranscripts,$n;
    push @altSquared,$n*$n; }
  if($code eq "mapped") {
    if($status->hasDescendent("NMD")) { ++$NMD{$gene} }
    if($status->hasDescendent("protein-differs")) { ++$proteinDiffers{$gene} }
    if($status->hasDescendent("frameshift")) {
      ++$frameshift{$gene};
      my $frameshift=$status->findDescendent("frameshift");
      my $length=0+$frameshift->getAttribute("nt-with-phase-mismatch");
      my $squaredLen=$length*$length;
      push @frameshifts,$length;
      push @squaredFrameshifts,$squaredLen;
    }
    if($status->hasDescendent("nonstop-decay")) { ++$nonstopDecay{$gene} }
    if($status->hasDescendent("premature-stop")) { ++$prematureStop{$gene} }
    if($status->hasDescendent("start-codon-change"))
      { ++$startCodonChange{$gene} }
  }
  elsif($code eq "no-transcript") {
    if($status->hasDescendent("broken-donor") ||
       $status->hasDescendent("broken-acceptor"))
      { ++$noTranscriptBrokenSS{$gene} }
    else { ++$noTranscriptWeakenedSS{$gene} }
  }
  elsif($code eq "splicing-changes" &&
	($status->hasDescendent("broken-donor") || 
	 $status->hasDescendent("broken-acceptor"))) {
    if($status->countDescendents("broken-donor")+
       $status->countDescendents("broken-acceptor")==1) {
      my $x=$status->countDescendentOrDatum("cryptic-site");
      $crypticSites+=$x; $crypticSitesSquared+=$x*$x;
      ++$singleBrokenSites;
      $x=$status->countDescendentOrDatum("exon-skipping");
      $exonSkipping+=$x; $exonSkippingSquared+=$x*$x;
    }
    if($status->hasDescendent("NMD")) { ++$splicingChangesNMD{$gene}; }
    if($coding) { ++$splicingChangesCoding{$gene} }
    else { ++$splicingChangesNoncoding{$gene} }
  }
}

my $NMD=keys %NMD; my $noTranscriptBrokenSS=keys %noTranscriptBrokenSS; 
my $noTranscriptWeakenedSS=keys %noTranscriptWeakenedSS;
my $splicingChangesNMD=keys %splicingChangesNMD; 
my $splicingChangesCoding=keys %splicingChangesCoding;
my $splicingChangesNoncoding=keys %splicingChangesNoncoding;
my $proteinDiffers=keys %proteinDiffers; my $frameshift=keys %frameshift;
my $nonstopDecay=keys %nonstopDecay; my $prematureStop=keys %prematureStop;
my $startCodonChange=keys %startCodonChange;
my $frameshiftCount=@frameshifts; my $totalFrameshiftLen=0;
foreach my $x (@frameshifts) {$totalFrameshiftLen+=$x}
my $totalFrameshiftSquared;
foreach my $x (@squaredFrameshifts) {$totalFrameshiftSquared+=$x}
my $numGenesWithAlt=@altTranscripts; my $totalAltTranscripts=0;
foreach my $x (@altTranscripts) {$totalAltTranscripts+=$x }
my $totalAltSquared;
foreach my $x (@altSquared) {$totalAltSquared+=$x}

if($VERBOSE) {
  print "#genes with any change to the protein: $proteinDiffers\n";
  print "#genes with a change in location of start codon: $startCodonChange\n";
  print "#transcripts that lost a start codon: $noStartCodon\n";
  print "#genes with a frameshift: $frameshift\n";
  print "#genes with nonstop decay (missing stop codon): $nonstopDecay\n";
  print "#genes with premature stop codon: $prematureStop\n";
  print "#genes with predicted NMD: $NMD\n";
  print "splicing changes:\n";
  print "\t#noncoding genes: $splicingChangesNoncoding\n";
  print "\t#coding genes: $splicingChangesCoding\n";
  print "\t\t#genes with NMD: $splicingChangesNMD\n";
  print "\t#broken donor sites: $brokenDonor\n";
  print "\t#broken acceptor sites: $brokenAcceptor\n";
  print "\t#cryptic sites near a broken site: $crypticSites\n";
  print "\t#weakened donor sites: $weakenedDonor\n";
  print "\t#weakened acceptor sites: $weakenedAcceptor\n";
  print "#genes with no transcript:\n\tdue to broken SS: $noTranscriptBrokenSS
\tdue to weakened SS: $noTranscriptWeakenedSS\n"
}
else {
  print "GENES_PROTEIN_CHANGE\t$proteinDiffers\n";
  print "GENES_CHANGE_START_CODON\t$startCodonChange\n";
  print "CODING_TRANSCRIPTS_NO_START_CODON\t$noStartCodon\n";
  print "GENES_FRAMESHIFT\t$frameshift\t$totalFrameshiftLen\t$totalFrameshiftSquared\t$frameshiftCount\n";
  print "GENES_NONSTOP_DECAY\t$nonstopDecay\n";
  print "GENES_PTC\t$prematureStop\n";
  print "MAPPED_GENES_NMD\t$NMD\n";
  print "NONCODING_GENES_SPLICING_CHANGE\t$splicingChangesNoncoding\n";
  print "CODING_GENES_SPLICING_CHANGE\t$splicingChangesCoding\n";
  print "GENES_SPLICING_CHANGE_NMD\t$splicingChangesNMD\n";
  print "BROKEN_DONORS\t$brokenDonor\n";
  print "BROKEN_ACCEPTORS\t$brokenAcceptor\n";
  print "CRYPTIC_SITES\t$crypticSites\t$crypticSitesSquared\t$singleBrokenSites\n";
  print "EXON_SKIPPING\t$exonSkipping\t$exonSkippingSquared\t$singleBrokenSites\n";
  print "ALT_TRANSCRIPTS\t$totalAltTranscripts\t$totalAltSquared\t$numGenesWithAlt\n";
  print "WEAKENED_DONORS\t$weakenedDonor\n";
  print "WEAKENED_ACCEPTORS\t$weakenedAcceptor\n";
  print "GENES_NO_TRANSCRIPT_BROKEN_SS\t$noTranscriptBrokenSS\n";
  print "GENES_NO_TRANSCRIPT_WEAK_SS\t$noTranscriptWeakenedSS\n"
}


