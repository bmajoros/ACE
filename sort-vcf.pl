#!/usr/bin/env perl
use strict;
use ProgramName;

# Process command line
my $name=ProgramName::get();
die "$name <in.vcf.gz> <out.vcf.gz> <temp-dir>\n" unless @ARGV==3;
my ($infile,$outfile,$tempDir)=@ARGV;

# Split out the header lines
$infile=~/([^\/]+)$/ || die $infile;
my $nopath=$1;
my $headerFile="$tempDir/$nopath.header";
my $bodyFile="$tempDir/$nopath.body";
open(IN,"cat $infile | gunzip |");
open(HEADER,">$headerFile") || die $headerFile;
open(BODY,">$bodyFile") || die $bodyFile;
while(<IN>) {
  if(/^#/) { print HEADER $_ }
  else { print BODY $_ }
}
close(BODY);
close(HEADER);
close(IN);

# Sort just the body
my $sortedFile="$tempDir/$nopath.sorted";
system("sort -k1,1d -k2,2n -T $tempDir $bodyFile > $sortedFile");

# Add in the header lines
system("cat $headerFile $sortedFile | bgzip > $outfile");

# Clean up
unlink($headerFile); unlink($bodyFile); unlink($sortedFile);


#cat /data/common/1000_genomes/VCF/20130502/GRCh38/chr10-hg38.vcf.gz | gunzip | sort -k1,1d -k2,2n -T /data/common/1000_genomes/VCF/20130502/GRCh38/tmp | bgzip > sorted/chr10-hg38.vcf.gz ; tabix sorted/chr10-hg38.vcf.gz
