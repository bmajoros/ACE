#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import ProgramName
from FastaReader import FastaReader
from FastaWriter import FastaWriter
import TempFilename

ACEPLUS=os.environ["ACEPLUS"]
ALPHABET=("A","C","G","T")

def emit(begin,end,examples,consensuses,category,OUT):
    nAlpha=len(ALPHABET)
    for example in examples:
        (defline,seq)=example
        #print("range",begin,end)
        consensus=seq[80:82]
        if(consensus not in consensuses): continue
        print(category,end="\t",file=OUT)
        for position in range(begin,end):
            nuc=seq[position]
            for i in range(nAlpha):
                letter=ALPHABET[i]
                count=1 if letter==nuc else 0
                print(count,end="",file=OUT)
                sep="\t" if i+1<nAlpha or position+1<end else "\n"
                print(sep,end="",file=OUT)

def printHeader(begin,end,OUT):
    print("category\t",end="",file=OUT)
    nAlpha=len(ALPHABET)
    for pos in range(begin,end):
        for i in range(nAlpha):
            nuc=ALPHABET[i]
            print(str(pos-begin)+nuc,end="",file=OUT)
            sep="\t" if i+1<nAlpha or pos+1<end else "\n"
            print(sep,end="",file=OUT)
            
#def emitPseudocounts(OUT):
#    print(1,end="",file=OUT)
#    for nuc in ALPHABET: print("\t1",end="",file=OUT)
#    print(file=OUT)
#    print(0,end="",file=OUT)
#    for nuc in ALPHABET: print("\t1",end="",file=OUT)
#    print(file=OUT)

def train(begin,end,positives,negatives,consensuses,featuresFile,betasFile):
    OUT=open(featuresFile,"wt")
    printHeader(begin,end,OUT)
    emit(begin,end,positives,consensuses,1,OUT)
    emit(begin,end,negatives,consensuses,0,OUT)
    OUT.close()
    os.system(ACEPLUS+"/logistic-regression.R "+featuresFile+" "+
              ALPHA+" "+betasFile) #+" 2> /dev/null")
    #exit(featuresFile+"\t"+betasFile)
    betas=loadBetas(betasFile)
    output(betas)

def output(betas):
    for pair in betas:
        (feature,beta)=pair
        print(feature,beta,sep="\t")

def loadBetas(filename):
    betas=[]
    with open(filename,"rt") as IN:
        IN.readline()
        IN.readline()
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): raise Exception(line)
            (feature,beta)=fields
            if(beta=="."): beta=0.0
            if(feature=="(Intercept)"): feature="intercept"
            else: feature=feature[1:]
            beta=round(float(beta),3)
            betas.append([feature,beta])
    return betas

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=7):
    exit(ProgramName.get()+" <pos.fasta> <neg.fasta> <left-context> <right-context> <alpha> <consensus,consensus,...\n")
(posFasta,negFasta,LEFT_MARGIN,RIGHT_MARGIN,ALPHA,consensuses)=sys.argv[1:]
consensuses=consensuses.split(",")
LEFT_MARGIN=int(LEFT_MARGIN)
RIGHT_MARGIN=int(RIGHT_MARGIN)

positives=FastaReader.readAllIntoArray(posFasta)
negatives=FastaReader.readAllIntoArray(negFasta)
tempFile1=TempFilename.generate(".fasta")
tempFile2=TempFilename.generate(".betas")
begin=80-LEFT_MARGIN
end=82+RIGHT_MARGIN
train(begin,end,positives,negatives,consensuses,tempFile1,tempFile2)
os.remove(tempFile1); os.remove(tempFile2)


