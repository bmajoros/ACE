#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
from GffTranscriptReader import GffTranscriptReader
from Translation import Translation
from Pipe import Pipe
from Interval import Interval
from FastaWriter import FastaWriter
import sys

def getSeq(chrom,strand,interval,twoBitFile):
    cmd="twoBitToFa -seq="+chrom+" -start="+str(interval.begin)+\
        " -end="+str(interval.end)+" "+twoBitFile+" stdout"
    pipe=Pipe(cmd)
    seq=""
    while(True):
        line=pipe.readline()
        if(not line): break
        if(len(line)>0 and line[0]=='>'): continue
        line=line.rstrip()
        seq+=line
    seq=seq.upper()
    if(strand=="-"): seq=Translation.reverseComplement(seq)
    return seq

# Process command line
if(len(sys.argv)!=5):
    exit(sys.argv[0]+" <in.gff> <genome.2bit> <out-exons.fasta> <out-introns.fasta")
(gffFile,twoBitFile,outExons,outIntrons)=sys.argv[1:]

# Load GFF
reader=GffTranscriptReader()
writer=FastaWriter()
EXONS=open(outExons,"wt")
INTRONS=open(outIntrons,"wt")
genes=reader.loadGenes(gffFile)
seen=set()
nextExonID=1
nextIntronID=1
for gene in genes:
    transcript=gene.longestTranscript()
    if(transcript.numExons()>0): continue # coding gene
    extra=transcript.parseExtraFields() # array of [key,value] pairs
    hash=transcript.hashExtraFields(extra)
    if(hash["gene_type"]!="lincRNA"): continue
    if(hash["gene_status"]!="KNOWN"): continue
    geneID=gene.getID()
    if(geneID in seen): continue
    seen.add(geneID)
    chrom=transcript.getSubstrate()
    strand=transcript.getStrand()
    for i in range(transcript.numUTR()):
        exon=transcript.getIthUTR(i)
        interval=Interval(exon.getBegin(),exon.getEnd())
        seq=getSeq(chrom,strand,interval,twoBitFile)
        if('N' in seq): continue
        defline=">"+str(nextExonID)+" /gene="+geneID+" /transcript="+\
            transcript.getID()+" /strand="+strand+" /exon="+str(i)
        nextExonID+=1
        writer.addToFasta(defline,seq,EXONS)
    introns=transcript.getIntrons()
    for intron in introns:
        seq=getSeq(chrom,strand,intron,twoBitFile)
        if('N' in seq): continue
        defline=">"+str(nextIntronID)+" /gene="+geneID+" /transcript="+\
            transcript.getID()+" /strand="+strand+" /intron="+str(i)
        nextIntronID+=1
        writer.addToFasta(defline,seq,INTRONS)
EXONS.close()
INTRONS.close()



