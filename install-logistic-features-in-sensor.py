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
from Translation import Translation
from NgramIterator import NgramIterator
from Rex import Rex
rex=Rex()

HALF_BETAS=True

def printHeader(featureType,N,OUT):
    print("IMM\n"+featureType+"\n"+str(N-1)+"\t-1\n"+str(N),file=OUT)

def writeOnes(order,OUT):
    print(str(4**(order+1)),file=OUT)
    ngramIterator=NgramIterator("ACGT",order+1)
    while(True):
        nmer=ngramIterator.nextString() # returns None if no more
        if(nmer is None): break
        print(nmer+"\n0.0",file=OUT) # 0.0 = log(1.0)

def writeModel(contentType,N,features,sign,revcomp,OUT):
    printHeader(contentType,N,OUT)
    for order in range(N-1): writeOnes(order,OUT)
    #for feature in features.keys():
    ngramIterator=NgramIterator("ACGT",N)
    print(str(4**N),file=OUT)
    while(True):
        feature=ngramIterator.nextString()
        if(feature is None): break
        score=sign*features.get(feature,0.0)
        if(revcomp): feature=Translation.reverseComplement(feature)
        print(feature+"\n"+str(score),file=OUT)

def loadFeatures(featureFile):
    global N
    features={}
    with open(featureFile,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (feature,score)=fields
            if(rex.find("Intercept",feature)): continue
            if(score=="."): score=0.0
            else: score=float(score)
            if(HALF_BETAS): score/=2.0
            features[feature]=score
            N=len(feature)
    return features

def normalize(features):
    maxValue=0.0
    for x in features.values():
        absX=abs(x)
        if(absX>maxValue): maxValue=absX
    for feature in features.keys():
        features[feature]/=maxValue

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+
         " <logistic-features.txt> <EXON|INTRON> <out.model>\n")
(featureFile,contentType,outFile)=sys.argv[1:]
contentTypeRev=None
if(contentType=="EXON"):
    contentType="SINGLE-EXON"
    contentTypeRev="NEG-SINGLE-EXON"
    sign=1
elif(contentType=="INTRON"):
    contentTypeRev="NEG-INTRON"
    sign=-1
else: raise Exception(contentType+" must be EXON or INTRON")

N=0
features=loadFeatures(featureFile)
#normalize(features)

OUT=open(outFile,"wt")
writeModel(contentType,N,features,sign,False,OUT)
writeModel(contentTypeRev,N,features,sign,True,OUT)
OUT.close()



