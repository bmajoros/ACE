/****************************************************************
 TataCapModel.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "TataCapModel.H"
#include <string.h>
#include <iostream>
#include <fstream>
#include "BOOM/Constants.H"
using namespace std;


TataCapModel::TataCapModel(TataCapModel &other,bool revComp)
  : SignalSensor(other.getGC(),other,revComp),
    minSeparation(other.minSeparation),
    maxSeparation(other.maxSeparation),
    tataLength(other.tataLength),
    capLength(other.capLength)
{
  if(revComp)
    {
      tataModel=(WMM*) other.tataModel->reverseComplement();
      capModel=(WMM*) other.capModel->reverseComplement();
      capIntergenicRatioModel=(WMM*)
	other.capIntergenicRatioModel->reverseComplement();
      intergenicModel=(MarkovChain*)
	other.intergenicModel->reverseComplement();
    }
  else
    {
      tataModel=new WMM(getGC(),*other.tataModel);
      capModel=new WMM(getGC(),*other.capModel);
      capIntergenicRatioModel=
	new WMM(getGC(),*other.capIntergenicRatioModel);
      intergenicModel=new MarkovChain(*other.intergenicModel);

      int contextWindowLength=tataLength+capLength+maxSeparation;
      int tataLen=tataModel->getConsensusLength();
      int tataOffset=tataModel->getConsensusOffset();
      setSizes(tataLen,
	       contextWindowLength-tataOffset-tataLen,
	       contextWindowLength);
    }
}



TataCapModel::TataCapModel(GarbageCollector &gc,
			   const BOOM::String &filename)
  : SignalSensor(gc)
{
  load(filename);
}



TataCapModel::TataCapModel(GarbageCollector &gc,istream &is)
  : SignalSensor(gc)
{
  load(is);
}



SignalSensor *TataCapModel::reverseComplement()
{
  return new TataCapModel(*this,true);
}



TataCapModel::~TataCapModel()
{
  delete tataModel;
  delete capModel;
  delete capIntergenicRatioModel;
  delete intergenicModel;
}



void TataCapModel::load(istream &is)
{
  is >> minSeparation >> maxSeparation;
  SignalType signalType=TSS;
  setSignalType(signalType);
  setStrand(::getStrand(signalType));
  tataModel=(WMM*) SignalSensor::load(is,getGC());
  intergenicModel=(MarkovChain*) ContentSensor::load(is);
  capModel=(WMM*) SignalSensor::load(is,getGC());
  capIntergenicRatioModel=(WMM*) SignalSensor::load(is,getGC());
  tataLength=tataModel->getContextWindowLength();
  capLength=capModel->getContextWindowLength();
  setSizes(tataModel->getConsensusLength(),
	   tataModel->getConsensusOffset(),
	   tataLength+capLength+maxSeparation);
  setCutoff(NEGATIVE_INFINITY); // ###
}



void TataCapModel::load(const BOOM::String &filename)
{
  ifstream is(filename.c_str());
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="TataCapModel")
    throw BOOM::String("Attempt to load a model of type ")+modelType
      +" into TataCapModel";
  load(is);
}



bool TataCapModel::consensusOccursAt(const BOOM::String &seq,int index)
{
  const char *str=seq.c_str();
  if(getStrand()==FORWARD_STRAND)
    {
      // First, there must at least be a TATA at the appropriate position in
      // the TATA window:
      if(strncmp(str+index,"TATA",4)) return false;

      // There must also be a C in the second position of the CAP window,
      // but since the CAP window can occur at a range of positions, we must
      // consider all of them:
      int tataWindowBegin=index-tataModel->getConsensusOffset();
      int spacerBegin=tataWindowBegin+tataLength;
      int firstCapWindowBegin=spacerBegin+minSeparation;
      int lastCapWindowBegin=spacerBegin+maxSeparation;
      for(int capWindowBegin=firstCapWindowBegin ; 
	  capWindowBegin<=lastCapWindowBegin ; 
	  ++capWindowBegin)
	{
	  int Cpos=capWindowBegin+1;
	  if(str[Cpos]=='C')
	    {
	      // Also, the base after the C must be an A or T:
	      char b=str[Cpos+1];
	      if(b=='A' || b=='T') 
		if(str[Cpos+3]!='G') 
		  return true;
	    }
	}
    }
  else
    {
      // First, there must at least be a TATA at the appropriate position in
      // the TATA window:
      index-=getConsensusOffset();
      int tataWindowBegin=index+capLength+maxSeparation;
      int consensusBegin=tataWindowBegin+tataLength-5;
      if(strncmp(str+consensusBegin,"TATA",4)) return false;

      // There must also be a C in the second position of the CAP window,
      // but since the CAP window can occur at a range of positions, we must
      // consider all of them:
      int lastCapWindowBegin=tataWindowBegin-minSeparation-capLength;
      for(int capWindowBegin=index ; 
	  capWindowBegin<=lastCapWindowBegin ; 
	  ++capWindowBegin)
	{
	  int lastCapPos=capWindowBegin+capLength-1;
	  int Cpos=lastCapPos-1;
	  if(str[Cpos]=='G')
	    {
	      char b=str[Cpos-1];
	      if(b=='A' || b=='T')
		if(str[Cpos-3]!='C') 
		  return true;
	    }
	}
    }
  return false;
}



double TataCapModel::getLogP(const Sequence &seq,const BOOM::String &str,
			     int begin)
{
  double logP;
  const char *cptr=str.c_str();
  if(getStrand()==FORWARD_STRAND)
    {
      // First, we need to evaluate the capIntergenicRatioModel at each 
      // possible position and choose the highest-scoring:
      int tataWindowBegin=begin;
      int spacerBegin=tataWindowBegin+tataLength;
      int firstCapWindowBegin=spacerBegin+minSeparation;
      int lastCapWindowBegin=spacerBegin+maxSeparation;
      int bestCapWindowBegin;
      double bestCapScore=NEGATIVE_INFINITY;;
      for(int capWindowBegin=firstCapWindowBegin ; 
	  capWindowBegin<=lastCapWindowBegin ; 
	  ++capWindowBegin)
	{
	  int Cpos=capWindowBegin+1;
	  if(cptr[Cpos]!='C') continue;
	  char b=cptr[Cpos+1];
	  if(b!='A' && b!='T' || cptr[Cpos+3]=='G') continue;
	  double capScore=
	    capIntergenicRatioModel->getLogP(seq,str,capWindowBegin);
	  if(capScore>bestCapScore || isinf(bestCapScore))
	    {
	      bestCapScore=capScore;
	      bestCapWindowBegin=capWindowBegin;
	    }
	}

      // Evaluate the TATA and CAP models at their respective positions:
      logP=tataModel->getLogP(seq,str,tataWindowBegin);
      logP+=capModel->getLogP(seq,str,bestCapWindowBegin);

      // Evaluate the intergenic spacer between the TATA & CAP models:
      for(int pos=spacerBegin ; pos<bestCapWindowBegin ; ++pos)
	logP+=
	  intergenicModel->scoreSingleBase(seq,str,pos,seq[pos],cptr[pos]);

      // Finally, evaluate the last few noncoding bases hanging on at
      // the end of the window after the CAP model (if any):
      int tailBegin=bestCapWindowBegin+capLength;
      int tailEnd=begin+getContextWindowLength();
      for(int pos=tailBegin ; pos<tailEnd ; ++pos)
	logP+=
	  intergenicModel->scoreSingleBase(seq,str,pos,seq[pos],cptr[pos]);
    }
  else
    {
      // First, we need to evaluate the capIntergenicRatioModel at each 
      // possible position and choose the highest-scoring:
      int tataWindowBegin=begin+capLength+maxSeparation;
      int lastCapWindowBegin=tataWindowBegin-minSeparation-capLength;
      int bestCapWindowBegin=-10000000;
      double bestCapScore=NEGATIVE_INFINITY;
      for(int capWindowBegin=begin ; 
	  capWindowBegin<=lastCapWindowBegin ; 
	  ++capWindowBegin)
	{
	  int lastCapPos=capWindowBegin+capLength-1;
	  int Cpos=lastCapPos-1;
	  if(cptr[Cpos]!='G') continue;
	  char b=cptr[Cpos-1];
	  if(b!='A' && b!='T' || cptr[Cpos-3]=='C') continue;
	  double capScore=
	    capIntergenicRatioModel->getLogP(seq,str,capWindowBegin);
	  if(capScore>bestCapScore || isinf(bestCapScore))
	    {
	      bestCapScore=capScore;
	      bestCapWindowBegin=capWindowBegin;
	    }
	}

      // Evaluate the TATA and CAP models at their respective positions:
      logP=tataModel->getLogP(seq,str,tataWindowBegin);
      logP+=capModel->getLogP(seq,str,bestCapWindowBegin);

      // Evaluate the intergenic spacer between the TATA & CAP models:
      int spacerBegin=bestCapWindowBegin+capLength;
      for(int pos=spacerBegin ; pos<tataWindowBegin ; ++pos)
	logP+=
	  intergenicModel->scoreSingleBase(seq,str,pos,seq[pos],cptr[pos]);

      // Finally, evaluate the last few noncoding bases hanging on at
      // the end of the window after the CAP model (if any):
      for(int pos=begin ; pos<bestCapWindowBegin ; ++pos)
	logP+=
	  intergenicModel->scoreSingleBase(seq,str,pos,seq[pos],cptr[pos]);
    }

  return logP;
}


