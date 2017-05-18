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

MAX_LEN=5000

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
if(len(sys.argv)!=7):
    exit(sys.argv[0]+" <in.gff> <genome.2bit> <out-exons.fasta> <out-introns.fasta> <skip#> <max-features>")
(gffFile,twoBitFile,outExons,outIntrons,skipNum,maxFeatures)=sys.argv[1:]
maxFeatures=int(maxFeatures); skipNum=int(skipNum)

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
    #if(transcript.numExons()>0): continue # skip coding genes
    #print(transcript.numExons(),"exons",flush=True)
    if(transcript.numExons()==0): continue # skip noncoding genes
    extra=transcript.parseExtraFields() # array of [key,value] pairs
    hash=transcript.hashExtraFields(extra)
    if(hash.get("gene_status") is not None and
       hash["gene_status"]!="KNOWN"): continue
    geneID=gene.getID()
    if(geneID in seen): continue
    seen.add(geneID)
    chrom=transcript.getSubstrate()
    strand=transcript.getStrand()
    #rawExons=transcript.getRawExons()
    rawExons=transcript.exons
    #print(len(transcript.UTR),"UTRs")
    for i in range(1,len(rawExons)-1):
        #print("i=",i,flush=True)
        if(nextExonID<skipNum): nextExonID+=1; continue
        if(nextExonID>maxFeatures): break
        exon=rawExons[i]
        interval=Interval(exon.getBegin(),exon.getEnd())
        seq=getSeq(chrom,strand,interval,twoBitFile)
        if('N' in seq): continue
        defline=">"+str(nextExonID)+" /gene="+geneID+" /transcript="+\
            transcript.getID()+" /strand="+strand+" /exon="+str(i)+\
            " /frame="+str(exon.getFrame())
        nextExonID+=1
        print("adding exon",len(seq),flush=True)
        writer.addToFasta(defline,seq,EXONS)
    introns=transcript.getIntrons()
    #print(len(introns),"introns",flush=True)
    for i in range(len(introns)):
        if(nextIntronID<skipNum): nextIntronID+=1; continue
        if(nextIntronID>maxFeatures): break
        intron=introns[i]
        #print("getSeq() ",intron.length(),intron.begin,intron.end,flush=True)
        seq=getSeq(chrom,strand,intron,twoBitFile)
        #print("checking Ns",flush=True)
        if('N' in seq): continue
        defline=">"+str(nextIntronID)+" /gene="+geneID+" /transcript="+\
            transcript.getID()+" /strand="+strand+" /intron="+str(i+1)
        nextIntronID+=1
        #print("adding to fasta",flush=True)
        writer.addToFasta(defline,seq,INTRONS)
        #print("done adding",flush=True)
    #print("end of loop")
    if(nextExonID>maxFeatures and nextIntronID>maxFeatures): break
EXONS.close()
INTRONS.close()



