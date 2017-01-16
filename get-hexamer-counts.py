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
import sys
import ProgramName
from FastaReader import FastaReader
from NgramIterator import NgramIterator

def getCounts(seq):
    hash={}
    for i in range(len(seq)-5):
        hexamer=seq[i:i+6]
        hash[hexamer]=hash.get(hexamer,0)+1
    #for hex in hash.keys(): print(hex,hash[hex])
    return hash

def getVector(counts,L):
    sampleSize=L-5
    vector=[]
    iter=NgramIterator("ACGT",6)
    while(True):
        hexamer=iter.nextString()
        if(hexamer is None): break
        count=counts.get(hexamer,0)
        fraction=float(count)/float(sampleSize)
        fraction=round(fraction,4)
        vector.append(fraction)
    return vector

def emit(vector,label):
    print(label,"\t",sep="",end="")
    for i in range(len(vector)):
        print(vector[i],"\t",sep="",end="")
    print("\n",end="")

def process(filename,label):
    reader=FastaReader(filename)
    while(True):
        (defline,seq)=reader.nextSequence()
        if(not defline): break
        if(len(seq)<50): continue
        counts=getCounts(seq)
        vector=getVector(counts,len(seq))
        emit(vector,label)
    reader.close()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <exons.fasta> <introns.fasta>\n")
(exonFile,intronFile)=sys.argv[1:]

# print header
print("category\t",end="")
iter=NgramIterator("ACGT",6)
while(True):
    hexamer=iter.nextString()
    if(hexamer is None): break
    print(hexamer,"\t",sep="",end="")
print("\n",end="")

process(exonFile,1)
process(intronFile,0)


