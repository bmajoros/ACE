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
import math
import sys
import os
import ProgramName
from FastaReader import FastaReader
from FastaWriter import FastaWriter
import TempFilename
from Rex import Rex
rex=Rex()

PSEUDOCOUNT=0.1
MIN_SENSITIVITY=0.99
WANT_TEST=False
POS_TEST_FILE="logreg-test.pos"
NEG_TEST_FILE="logreg-test.neg"
POS_RAW_FILE="logreg-raw.pos"
NEG_RAW_FILE="logreg-raw.neg"
TEST_ROC="logreg-test.roc"
ACEPLUS=os.environ["ACEPLUS"]
ALPHABET=("A","C","G","T")
QUIET=True

def getModel(betas):
    model={}
    for elem in betas:
        (feature,beta)=elem
        model[feature]=beta
    return model

def getP(model,seq):
    score=getRawScore(model,seq)
    P=1.0/(1.0+math.exp(-score))
    return P

def getRawScore(model,seq):
    L=len(seq)
    score=model.get("intercept",0.0)
    for i in range(L):
        key=str(i)+seq[i]
        beta=model[key]
        score+=beta
    return score

def getScores(model,examples,begin,end):
    scores=[]
    for example in examples:
        (defline,seq)=example
        subseq=seq[begin:end]
        P=getP(model,subseq)
        scores.append(P)
    return scores

def getRawScores(model,examples,begin,end):
    scores=[]
    for example in examples:
        (defline,seq)=example
        subseq=seq[begin:end]
        score=getRawScore(model,subseq)
        scores.append(score)
    return scores

def test(betas,positives,negatives,begin,end):
    model=getModel(betas)
    posScores=getScores(model,positives,begin,end)
    negScores=getScores(model,negatives,begin,end)
    writeScores(posScores,POS_TEST_FILE)
    writeScores(negScores,NEG_TEST_FILE)
    OUT=open(TEST_ROC,"wt")
    writeROC(posScores,1,OUT)
    writeROC(negScores,0,OUT)
    OUT.close()
    posScores=getRawScores(model,positives,begin,end)
    negScores=getRawScores(model,negatives,begin,end)
    writeScores(posScores,POS_RAW_FILE)
    writeScores(negScores,NEG_RAW_FILE)

def writeROC(scores,label,OUT):
    for score in scores:
        print(score,label,sep="\t",file=OUT)

def writeScores(scores,filename):
    with open(filename,"wt") as OUT:
        for score in scores:
            print(score,file=OUT)

def emit(begin,end,examples,consensuses,category,OUT):
    nAlpha=len(ALPHABET)
    for example in examples:
        (defline,seq)=example
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

def getFreqs(begin,end,examples,consensuses):
    counts={}
    nAlpha=len(ALPHABET)
    for example in examples:
        (defline,seq)=example
        consensus=seq[80:82]
        if(consensus not in consensuses): continue
        for position in range(begin,end):
            nuc=seq[position]
            key=str(position-begin)+nuc
            counts[key]=counts.get(key,0)+1
    for pos in range(0,end-begin):
        sum=0.0
        for nuc in ALPHABET:
            if(nuc=="N"): continue
            key=str(pos)+nuc
            sum+=float(counts.get(key,PSEUDOCOUNT))
    for pos in range(0,end-begin):
        for nuc in ALPHABET:
            if(nuc=="N"): continue
            key=str(pos)+nuc
            f=float(counts.get(key,PSEUDOCOUNT))/sum
            counts[key]=math.log(f)
    return counts

def printHeader(begin,end,OUT):
    print("category\t",end="",file=OUT)
    nAlpha=len(ALPHABET)
    for pos in range(begin,end):
        for i in range(nAlpha):
            nuc=ALPHABET[i]
            print(str(pos-begin)+nuc,end="",file=OUT)
            sep="\t" if i+1<nAlpha or pos+1<end else "\n"
            print(sep,end="",file=OUT)
            
def train(begin,end,positives,negatives,consensuses,featuresFile,betasFile):
    OUT=open(featuresFile,"wt")
    printHeader(begin,end,OUT)
    emit(begin,end,positives,consensuses,1,OUT)
    emit(begin,end,negatives,consensuses,0,OUT)
    OUT.close()
    cmd=ACEPLUS+"/logistic-regression.R "+featuresFile+" "+\
        ALPHA+" "+betasFile
    if(QUIET): cmd+=" 2> /dev/null"
    os.system(cmd)
    betas=loadBetas(betasFile)
    return betas

def getCutoff(betas,positives,begin,end):
    model=getModel(betas)
    scores=getRawScores(model,positives,begin,end)
    scores.sort()
    threshold=round(scores[int(len(scores)*(1.0-MIN_SENSITIVITY))],3)
    return threshold

def getFreqCutoff(freqs,positives,begin,end):
    model=freqs
    scores=getRawScores(model,positives,begin,end)
    scores.sort()
    threshold=round(scores[int(len(scores)*(1.0-MIN_SENSITIVITY))],3)
    return threshold

def output(betas,signalType,consensusOffset,consensusLen,windowLen,threshold,
           freqs):
    print("LogisticSensor")
    alphabetSize=5 # ACGNT
    #threshold=0.0
    #for pair in betas:
    #    (feature,beta)=pair
    #    beta=float(beta)/2.0
    #    if(feature=="intercept"): 
    #        threshold=beta
    #        break
    print(signalType,threshold,windowLen,alphabetSize,windowLen,
          consensusOffset,consensusLen,"+",sep="\t")
    numPairs=len(betas)
    print(numPairs)
    for pair in betas:
        (feature,beta)=pair
        beta=float(beta)/2.0
        pos=0; nuc=""
        if(feature=="intercept"): 
            print(feature,beta,sep="\t")
        else:
            if(not rex.find("(\d+)(\S+)",feature)):
                raise Exception("can't parse "+feature)
            print(rex[1],rex[2],beta,sep="\t")
    print()
    #for pos in range(windowLen):
    #    for nuc in ALPHABET:
    #        if(nuc=="N"): continue
    #        key=str(pos)+nuc
    #        value=freqs.get(key,0.0)
    #        value=round(value,4)
    #        print(pos,nuc,value,sep="\t")

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
if(len(sys.argv)!=8):
    exit(ProgramName.get()+" <pos.fasta> <neg.fasta> <left-context> <right-context> <alpha> <consensus-list> <type=GT|AG>\n")
(posFasta,negFasta,LEFT_MARGIN,RIGHT_MARGIN,ALPHA,consensuses,signalType)=\
    sys.argv[1:]
consensuses=consensuses.split(",")
LEFT_MARGIN=int(LEFT_MARGIN)
RIGHT_MARGIN=int(RIGHT_MARGIN)

positives=FastaReader.readAllIntoArray(posFasta)
negatives=FastaReader.readAllIntoArray(negFasta)
tempFile1=TempFilename.generate(".fasta")
tempFile2=TempFilename.generate(".betas")
begin=80-LEFT_MARGIN
end=82+RIGHT_MARGIN
betas=train(begin,end,positives,negatives,consensuses,tempFile1,tempFile2)
os.remove(tempFile1); os.remove(tempFile2)
consensusOffset=LEFT_MARGIN
consensusLen=2 # splice sites only!
windowLen=LEFT_MARGIN+2+RIGHT_MARGIN
freqs=getFreqs(begin,end,positives,consensuses)
cutoff=getCutoff(betas,positives,begin,end)
#cutoff=getFreqCutoff(freqs,positives,begin,end)
output(betas,signalType,consensusOffset,consensusLen,windowLen,cutoff,freqs)
if(WANT_TEST): test(betas,positives,negatives,begin,end)

