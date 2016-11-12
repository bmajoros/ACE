/****************************************************************
 FastMarkovChain.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "FastMarkovChain.H"
#include <math.h>
#include <iostream>
#include <fstream>
#include "BOOM/ProteinTrans.H"
#include "MarkovChainCompiler.H"

extern Alphabet alphabet;



FastMarkovChain::FastMarkovChain(BOOM::File &file)
  : alphabetSize(alphabet.getNumElements()),
    state(0),
    transitionMatrix(0,0),
    probabilityVector(0),
    statesPerOrder(0),
    revComp(NULL)
{
  load(file);
}



FastMarkovChain::FastMarkovChain(const BOOM::String &filename)
  : alphabetSize(alphabet.getNumElements()),
    state(0),
    transitionMatrix(0,0),
    probabilityVector(0),
    statesPerOrder(0),
    revComp(NULL)
{
  load(filename);
}



FastMarkovChain::FastMarkovChain(int order,ContentType contentType)
  : transitionMatrix(0,0),
    probabilityVector(0),
    statesPerOrder(order+1),
    order(order),
    alphabetSize(alphabet.getNumElements()),
    state(0),
    revComp(NULL)
{
  setContentType(contentType);
  computeNumberOfStates();

  probabilityVector.resize(numStates);
  transitionMatrix.resize(numStates,alphabetSize);  
}



void FastMarkovChain::computeNumberOfStates()
{
  numStates=1; // initial state corresponds with no combination
  int numCombinations=1; // allow N^order combinations
  for(int i=0 ;i<=order ; ++i)
    {
      numCombinations*=alphabetSize; // N^order
      statesPerOrder[i]=numCombinations;
      numStates+=numCombinations;
    }
}



bool FastMarkovChain::save(const BOOM::String &filename)
{
  BOOM::File file(filename,"w");
  return save(file);
}



void FastMarkovChain::reset(const Sequence &seq,const BOOM::String &str,
			    int pos)
{
  switch(getStrand())
    {
    case NO_STRAND:
    case FORWARD_STRAND:
      state=0;
      break;
    case REVERSE_STRAND:
      {
	seqLen=str.length();
	int len=order+1;
	if(pos+len>=seqLen) len=seqLen-pos;
	state=stringToState(str.substr(pos,len));
      }
      break;
    }
}



void FastMarkovChain::load(const BOOM::String &filename)
{
  BOOM::File file(filename);
  load(file);
}



void FastMarkovChain::load(BOOM::File &file)
{
  loadNonRecurs(file);
  ContentType contentType=getContentType();
  if(contentType!=INTERGENIC)
    {
      revComp=new FastMarkovChain(order,::reverseComplement(contentType));
      revComp->revComp=this;
      revComp->loadNonRecurs(file);
    }
}



void FastMarkovChain::loadNonRecurs(BOOM::File &file)
{
  BOOM::String modelType;
  file >> modelType;
  if(modelType!="FastMC")
    throw BOOM::String("Attempt to load model of type \"")+modelType+
      "\" into FastMarkovChain from binary file";

  int contentTypeCode;
  file >> order >> contentTypeCode;
  setContentType(static_cast<ContentType>(contentTypeCode));
  setStrand(::getStrand(static_cast<ContentType>(contentTypeCode)));
  statesPerOrder.resize(order+1);
  computeNumberOfStates();
  transitionMatrix.loadBytes(file);
  probabilityVector.loadBytes(file);
}



bool FastMarkovChain::save(BOOM::File &file)
{
  if(!saveNonRecurs(file)) return false;
  if(getContentType()!=INTERGENIC)
    {
      reverseComplement();
      return revComp->saveNonRecurs(file);
    }
  return true;
}



bool FastMarkovChain::saveNonRecurs(BOOM::File &file)
{
  file << "FastMC" << order << static_cast<int>(getContentType());
  transitionMatrix.saveBytes(file);
  probabilityVector.saveBytes(file);
  return true;
}



double FastMarkovChain::scoreSingleBase(const Sequence &seq,
					const BOOM::String &str,
					int index,Symbol s,char)
{
  switch(getStrand())
    {
    case NO_STRAND:
    case FORWARD_STRAND:
      state=transitionMatrix[state][s];
      break;
    case REVERSE_STRAND:
      {
	int nextIndex=index+order;
	if(nextIndex<seqLen)
	  {
	    state=transitionMatrix[state][seq[nextIndex]];
	  }
	else
	  {
	    state=stringToState(str.substr(index,seqLen-index));
	  }
      }
      break;
    }
  return probabilityVector[state];
}



void FastMarkovChain::scoreSingleBase(const Sequence &ss,
				      const BOOM::String &st,
				      int index,Symbol s,char c,
				      double &scorePhase0,
				      double &scorePhase1,
				      double &scorePhase2)
{
  scorePhase0=scorePhase1=scorePhase2=
    scoreSingleBase(ss,st,index,s,c);
}



double FastMarkovChain::scoreSubsequence(const Sequence &seq,
					 const BOOM::String &str,
					 int begin,int length,int phase)
{
  char c;
  double score=0;
  int end=begin+length;
  for(int pos=begin ; pos<end ; ++pos)
    score+=scoreSingleBase(seq,str,pos,seq[pos],c);
  return score;
}



bool FastMarkovChain::save(ostream &os)
{
  throw "attempt to save FastMarkovChain in non-binary file";
}



void FastMarkovChain::setTrans(int state1,Symbol s,int state2)
{
  transitionMatrix[state1][s]=state2;
}



void FastMarkovChain::setProb(int state,double prob)
{
  probabilityVector[state]=prob;
}



int FastMarkovChain::statesOfOrder(int order)
{
  return statesPerOrder[order];
}



int FastMarkovChain::stringToState(const BOOM::String &s)
{
  int len=s.length();
  if(len==0) return 0;
  int state=1;
  for(int i=0 ; i<len-1 ; ++i)
    state+=statesOfOrder(i);
  int alphabetSize=alphabet.getNumElements();
  int factor=1;
  for(int i=0 ; i<len ; ++i)
    {
      state+=(alphabet.lookup(s[len-i-1]) * factor);
      factor*=alphabetSize;
    }
  return state;
}



BOOM::String FastMarkovChain::stateToString(int state)
{
  int alphabetSize=alphabet.getNumElements();
  int order=0;
  int firstStateOfOrder=1;
  for(; 1 ; ++order)
    {
      if(state<firstStateOfOrder) break;
      firstStateOfOrder+=statesOfOrder(order);
    }
  --order;
  firstStateOfOrder=1;
  for(int i=0 ; i<order ; ++i)
    firstStateOfOrder+=statesOfOrder(i);
  state-=firstStateOfOrder;
  BOOM::String str;
  for(int i=0 ; i<=order ; ++i)
    {
      int factor=static_cast<int>(pow((float)alphabetSize,order-i));
      for(int s=alphabetSize-1 ; s>=0 ; --s)
	{
	  if(state>=s*factor)
	    {
	      str+=alphabet.lookup(s);
	      state-=s*factor;
	      break;
	    }
	}
    }
  return str;
}



BOOM::String FastMarkovChain::getState()
{
  return stateToString(state)+"("+state+")";
}



void FastMarkovChain::gotoState(const BOOM::String &stateLabel)
{
  state=stringToState(stateLabel);
}



ContentSensor *FastMarkovChain::reverseComplement()
{
  if(!revComp && getContentType()==INTERGENIC) revComp=this;
  return revComp;
}


