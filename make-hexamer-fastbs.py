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
import math

#SCALE=math.log(10)

def normalizeScores(scores):
    maxValue=0.0
    for pair in scores:
        (hexamer,score)=pair
        pair[1]=score=float(score)
        absValue=abs(score)
        if(absValue>maxValue): maxValue=absValue
    #maxValue/=SCALE
    for pair in scores:
        pair[1]/=maxValue

def writeScores(scores,category,filename):
    OUT=open(filename,"wt")
    for pair in scores:
        (hexamer,score)=pair
        score=math.exp(float(score)*category)
        print(hexamer,score,sep="\t",file=OUT)
    OUT.close()

#=========================================================================
# main()
#=========================================================================

if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <in:betas.txt> <out:pos-weights.txt> <out:neg-weights.txt> <outdir>\n")
(betasFile,outPos,outNeg,outDir)=sys.argv[1:]

scores=[]
with open(betasFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (hexamer,score)=fields
        if(hexamer=="(Intercept)"): continue
        if(score=="."): score=0
        scores.append([hexamer,score])
        outfile=outDir+"/"+hexamer+".fastb"
        OUT=open(outfile,"wt")
        OUT.write(">dna\n"+hexamer+"\n")
        OUT.close()
normalizeScores(scores)
writeScores(scores,1.0,outPos)
writeScores(scores,-1.0,outNeg)



