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
from Translation import Translation
from Rex import Rex
rex=Rex()

def readModels(fh,models):
    for line in fh:
        if(not line): break
        if(not rex.find("\S",line)): break
        fields=line.split()
        if(len(fields)!=2): exit("syntax error in model file: "+line)
        (kmer,score)=fields
        order=len(kmer)-1
        models[order][kmer]=score

def writeModel(hash,OUT):
    keys=hash.keys()
    N=len(keys)
    OUT.write(str(N)+"\n")
    for key in keys:
        score=hash[key]
        OUT.write(key+"\n")
        OUT.write(str(score)+"\n")

def writeModelReverse(hash,OUT):
    keys=hash.keys()
    N=len(keys)
    OUT.write(str(N)+"\n")
    for key in keys:
        score=hash[key]
        OUT.write(Translation.reverseComplement(key)+"\n")
        OUT.write(str(score)+"\n")

def writeModelFile(models,filename,contentType,order):
    with open(filename,"wt") as OUT:
        OUT.write("IMM\n")
        OUT.write(contentType+"\n")
        OUT.write(str(order)+"\t-1\n")
        OUT.write(str(order+1)+"\n")
        for i in range(order+1):
            writeModel(models[i],OUT)
        # also write reverse-strand model:
        OUT.write("IMM\n")
        if(contentType=="SINGLE-EXON"): contentType="NEG-SINGLE-EXON"
        elif(contentType=="INTRON"):  contentType="NEG-INTRON"
        OUT.write(contentType+"\n")
        OUT.write(str(order)+"\t-1\n")
        OUT.write(str(order+1)+"\n")
        for i in range(order+1):
            writeModelReverse(models[i],OUT)

#=========================================================================
#                                 main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(sys.argv[0]+" <in.mummie> <exon|intron|intergenic> <out.ace>")
(infile,contentType,outfile)=sys.argv[1:]
if(contentType=="exon"): contentType="SINGLE-EXON"
elif(contentType=="intergenic"): contentType="INTERGENIC"
elif(contentType=="intron"): contentType="INTRON"
else: exit("unknown content type: "+contentType)
order=None
models=[] # indexed by order
with open(infile,"rt") as IN:
    for line in IN:
        if(rex.find("(\d+) order",line)):
            order=int(rex[1])
            for i in range(order+1): models.append({})
            IN.readline() # alphabet:
            IN.readline() # ACGT
            readModels(IN,models)
            break
writeModelFile(models,outfile,contentType,order)

