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
import sys

# Process command line
if(len(sys.argv)!=2):
    exit(sys.argv[0]+" <in.gff>")
(gffFile,)=sys.argv[1:]

# Load GFF
reader=GffTranscriptReader()
genes=reader.loadGenes(gffFile)
sum=0
N=0
for gene in genes:
    transcript=gene.longestTranscript()
    extra=transcript.parseExtraFields() # array of [key,value] pairs
    hash=transcript.hashExtraFields(extra)
    if(hash.get("gene_status") is not None and
       hash["gene_status"]!="KNOWN"): continue
    rawExons=transcript.getRawExons()
    for i in range(1,len(rawExons)-1): # interal exons only
        exon=rawExons[i]
        sum+=exon.getLength()
        N+=1
mean=round(float(sum)/float(N))
print(mean)



