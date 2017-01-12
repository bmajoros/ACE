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
from FastaReader import FastaReader
from Rex import Rex
rex=Rex()

def computeScore(seq,hash):
    L=len(seq)
    end=L-6
    sum=0
    for i in range(end):
        hex=seq[i:i+6]
        sum+=hash[hex]
    score=sum/end
    return score

def loadHexamers(filename):
    hash={}
    with open(filename,"rt") as fh:
        for line in fh:
            line=line.rstrip()
            fields=line.split()
            (hex,score)=fields
            hash[hex]=float(score)
    return hash

if(len(sys.argv)!=4):
    exit(sys.argv[0]+" <in.fasta> <in.hexamers> <out.txt>")
(fastaFile,hexFile,outFile)=sys.argv[1:]

hash=loadHexamers(hexFile)

OUT=open(outFile,"wt")
reader=FastaReader(fastaFile)
while(True):
    (defline,seq)=reader.nextSequence()
    if(not defline): break
    score=computeScore(seq,hash)
    if(not rex.find("^\s*>\s*(\S+)",defline)):
        exit("can't parse defline: "+defline)
    id=rex[1]
    OUT.write(id+"\t"+str(score)+"\n")
OUT.close()    



