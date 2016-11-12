#!/bin/env perl
use strict;
use FileHandle;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;
use ConfigFile;
use Getopt::Std;
our ($opt_d);
getopts('d');

die "\n
major-allele-genome.pl [-d] <ace.config> <genes.gff> <out-dir> <population-list.txt> <indiv-population.txt>

-d = dry run: no output, just report errors

Assumptions:
 * populations.txt has two columns: individuals and populations
 * VCF files must be zipped with bgzip and have accompanying tbi indexes
 * VCF files must contain chromosome in filename, e.g. chr14, chrX, etc.
 * Chrom names begin with \"chr\" in all VCF and GF files
 * Your environment variable \$ACE must point to the ACE dir
 * <out-dir> will be populated with FASTA files
"
  unless @ARGV==5;
my ($configFile,$gffFile,$outDir,$IDfile,$populationFile)=@ARGV;

#==============================================================
# First, some initialization
#==============================================================


my $DRY_RUN=$opt_d;
my $DEBUG=0;
my $VERBOSE=1;
my $MARGIN_AROUND_GENE=1000;
my $config=new ConfigFile($configFile);
my $CHROM_LENGTHS=$config->lookupOrDie("chr-lengths");
my $TABIX=$config->lookupOrDie("tabix");
my $twoBitFile=$config->lookupOrDie("genome");
my $twoBitToFa=$config->lookupOrDie("twoBitToFa");
my $vcfDir=$config->lookupOrDie("vcf");
my $ploidy=0+$config->lookupOrDie("ploidy");
my $vcfLacksChr=$config->lookup("vcf-lacks-chr");
$vcfLacksChr=($vcfLacksChr eq "true");
if($ploidy<1) { die "invalid ploidy" }
my %chromLen;
loadChromLengths($CHROM_LENGTHS);
my %chrToVCF;
initChrToVCF($vcfDir);
system("mkdir -p $outDir") unless -e $outDir;
system("rm -f $outDir/errors.txt");
my $refGeneFasta="$outDir/refgene.fasta";
my $altGeneFasta="$outDir/altgene.fasta";
my $tempBedFile="$outDir/temp.bed";
my $geneVcfFile="$outDir/gene.vcf";#"$outDir/gene.vcf.gz";
my $geneTvfFile="$outDir/gene.tvf";#"$outDir/gene.tvf.gz";
my $populationVCF="$outDir/population.vcf";
my $outGFF="$outDir/local.gff";
my $ACE=$ENV{"ACE"};
my $fastaWriter=new FastaWriter;

#==============================================================
# Load gene coordinates from GFF file
#==============================================================

my $gffReader=new GffTranscriptReader();
my $genes=$gffReader->loadGenes($gffFile);

#==============================================================
# Make FASTA files for each individual
#==============================================================

my %keepIDs;
loadIDs($IDfile,\%keepIDs);
$keepIDs{"reference"}=1;
my @individuals=keys %keepIDs;
my $numIndiv=@individuals;
my %fastaFiles;
for(my $i=0 ; $i<$numIndiv ; ++$i) {
  my $indiv=$individuals[$i];
  next unless $keepIDs{$indiv};
  $fastaFiles{$indiv}=[];
  for(my $j=1 ; $j<=$ploidy ; ++$j) {
    my $file="$outDir/$indiv-$j.fasta";
    push @{$fastaFiles{$indiv}},$file;
  }
}
$fastaFiles{"reference"}=[];
for(my $j=1 ; $j<=$ploidy ; ++$j) {
  my $file="$outDir/ref-$j.fasta";
  push @{$fastaFiles{"reference"}},$file;
}

#==============================================================
# Process each gene
#==============================================================

my $numGenes=@$genes;
print "$numGenes genes loaded\n";
my %skipped;
if(!$DRY_RUN) { open(GFF,">$outGFF") || die "Can't create file $outGFF" }
for(my $i=0 ; $i<$numGenes ; ++$i) {
  my $gene=$genes->[$i];
  my $geneName=$gene->getId();
  print "gene $geneName\n";
  my $chr=$gene->getSubstrate();
  my $chrVcfFile=$chrToVCF{$chr};
  if(!$chrVcfFile) { $skipped{$chr}=1; print "Warning: skipping $chr (no VCF)\n"; next }
  my $begin=$gene->getBegin()-$MARGIN_AROUND_GENE;
  my $end=$gene->getEnd()+$MARGIN_AROUND_GENE;
  if(!$DRY_RUN) { writeLocalGFF($gene,$begin,*GFF) }
  if($begin<0) { $begin=0 }
  if(!defined($chromLen{$chr})) { $skipped{$chr}=1; print "Warning: skipping $chr (no length)\n"; next }
  if($end>$chromLen{$chr}) { $end=$chromLen{$chr} }
  my $strand=$gene->getStrand();
  my $name=$gene->getId();
  writeBed4($chr,$begin,$end,$name,$tempBedFile);
  System("$twoBitToFa -bed=$tempBedFile -noMask $twoBitFile $refGeneFasta");
  writeBed3($chr,$begin,$end,$tempBedFile);
  System("$TABIX -h $chrVcfFile -R $tempBedFile > $geneVcfFile");
  System("$ACE/vcf-population $geneVcfFile $populationFile $populationVCF");
  my $dashC=$vcfLacksChr ? " -c " : "";
  System("$ACE/vcf-to-tvf $dashC -i $IDfile $populationVCF $geneTvfFile");
  writeBed6($chr,$begin,$end,$name,$strand,$tempBedFile);
  system("rm -f $altGeneFasta");
  my $dashD=$DRY_RUN ? "-d" : "";
  my $errFile="$outDir/err.out";
  System("$ACE/tvf-to-fasta -p 1 $dashD -r $geneTvfFile $twoBitFile $tempBedFile $altGeneFasta >& $errFile");
  my $err=`cat $errFile`;
  if($err=~/error/ || $err=~/Abort/) { die $err }
  my (%warnings,%errors);
  loadErrors($errFile,\%warnings,\%errors);
  system("cat $errFile >> $outDir/errors.txt");
  next if $DRY_RUN;
  die unless -e $altGeneFasta;
  die if -z $altGeneFasta;
  my $fastaReader=new FastaReader($altGeneFasta);
  while(1) {
    my ($def,$seq)=$fastaReader->nextSequence();
    last unless $def;
    $def=~/>\S+\s+\/individual=(\S+)\s+\/allele=(\d+)\s+\/locus=(\S+)\s+\/coord=(\S+)\s+\/cigar=(\S+)\s+\/variants=(\S*)/
      || die "Can't parse defline: $def\n";
    my ($indivID,$alleleNum,$geneID,$coord,$cigar,$variants)=
      ($1,$2,$3,$4,$5,$6);
    if($keepIDs{$indivID}) {
      my $file=$fastaFiles{$indivID}->[$alleleNum-1];
      my $key="$indivID $geneID";
      my $numWarn=0+$warnings{$key};
      my $numErr=0+$errors{$key};
      open(FASTA,">>$file") || die $file;
      $def=">${geneID}_$alleleNum /coord=$coord /margin=$MARGIN_AROUND_GENE /cigar=$cigar /warnings=$numWarn /errors=$numErr /variants=$variants";
      $fastaWriter->addToFasta($def,$seq,\*FASTA);
      close(FASTA);
    }
    undef $seq; undef $def;
    undef $indivID; undef $alleleNum ; undef $geneID ; undef $coord;
  }
  $fastaReader->close();
  last if $DEBUG;
}
close(GFF);
my @skipped=keys %skipped;
foreach my $skipped (@skipped) { print "warning: skipped $skipped\n" }

#==============================================================
# Clean up
#==============================================================

if(!$DEBUG) {
  unlink($refGeneFasta);
  unlink($altGeneFasta);
  unlink($tempBedFile);
  unlink($geneVcfFile);
  unlink($geneTvfFile);
}

#==============================================================
# writeBed4($chr,$begin,$end,$name,$tempBedFile);
sub writeBed4 {
  my ($chr,$begin,$end,$name,$outfile)=@_;
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\t$name\n";
  close(OUT);
}
#==============================================================
# writeBed3($chr,$begin,$end,$tempBedFile);
sub writeBed3 {
  my ($chr,$begin,$end,$outfile)=@_;

  if($vcfLacksChr && $chr=~/chr(.+)/) { $chr=$1 }
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\n";
  close(OUT);
}
#==============================================================
# writeBed6($chr,$begin,$end,$name,$strand,$tempBedFile);
sub writeBed6 {
  my ($chr,$begin,$end,$name,$strand,$outfile)=@_;
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\t$name\t0\t$strand\n";
  close(OUT);
}
#==============================================================
sub System {
  my ($cmd)=@_;
  if($VERBOSE) { print "$cmd\n" }
  system($cmd);
}
#==============================================================
# initChrToVCF($vcfDir);
sub initChrToVCF {
  my ($dir)=@_;
  my @files=`ls $dir/*.vcf.gz`;
  foreach my $file (@files) {
    chomp $file;
    if($file=~/(chr[A-Za-z\d]+)/) {
      my $chr=$1;
      $chrToVCF{$chr}=$file;
    }
  }
}
#==============================================================
sub loadChromLengths
{
  my ($infile)=@_;
  open(IN,$infile) || die $infile;
  #<IN>; no header!
  while(<IN>) {
    chomp;
    my @fields=split;
    next unless @fields>=2;
    my ($chr,$len)=@fields;
    $chromLen{$chr}=$len;
  }
  close(IN);
}
#==============================================================
sub loadIDs
{
  my ($file,$hash)=@_;
  open(IN,$file) || die "can't open file $file\n";
  while(<IN>) {
    chomp;
    next unless $_=~/$\s*(\S+)/;
    my $ID=$1;
    $hash->{$ID}=1;
  }
  close(IN);
}
#==============================================================
sub loadErrors
{
  my ($filename,$warnings,$errors)=@_;
  open(IN,$filename) || die "can't open file: $filename";
  while(<IN>) {
    chomp;
    my @fields=split;
    next unless @fields>=6;
    my ($severity,$type,$indiv,$gene)=@fields;
    my $key="$indiv $gene";
    if($severity=~/WARNING/) { ++$warnings->{$key} }
    elsif($severity=~/ERROR/) { ++$errors->{$key} }
  }
  close(IN);
}
#==============================================================
sub writeLocalGFF
{
  my ($gene,$begin,$fh)=@_;
  my $numTranscripts=$gene->getNumTranscripts();
  for(my $i=0 ; $i<$numTranscripts ; ++$i) {
    my $transcript=$gene->getIthTranscript($i);
    $transcript->shiftCoords(-$begin);
    my $gff=$transcript->toGff();
    print $fh $gff;
  }
}
#==============================================================
#==============================================================

