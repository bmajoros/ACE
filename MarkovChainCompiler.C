/****************************************************************
 MarkovChainCompiler.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "MarkovChainCompiler.H"
#include <math.h>
#include <iostream>
#include "BOOM/Alphabet.H"
#include "NthOrderStringIterator.H"

extern Alphabet alphabet;
const double LOG_OF_1=0.0;



FastMarkovChain *MarkovChainCompiler::compile(ContentSensor &mc)
{
  FastMarkovChain *fmc=compileForward(mc);
  if(mc.getContentType()!=INTERGENIC)
    fmc->revComp=static_cast<FastMarkovChain*>
      (compileReverse(*mc.reverseComplement()));
  return fmc;
}



FastMarkovChain *MarkovChainCompiler::compileForward(ContentSensor &mc)
{
  int order=mc.getOrder();
  FastMarkovChain *fmc=new FastMarkovChain(order,mc.getContentType());
  int alphabetSize=alphabet.getNumElements();

  // ===========================
  // Install state probabilities
  // ===========================
  Sequence seq;
  Symbol s;
  fmc->setProb(0,LOG_OF_1);
  int state=1;
  for(int i=0 ; i<=order ; ++i)
    {
      NthOrderStringIterator stateGen(i+1,alphabet);
      while(!stateGen.done())
	{
	  BOOM::String stateString=stateGen.getNextString();
	  char c=stateString[i];
	  double score=mc.scoreSingleBase(seq,stateString,i,s,c);
	  fmc->setProb(state,score);
	  ++state;
	}
    }  

  // ===================
  // Install transitions
  // ===================

  for(int i=0 ; i<alphabetSize ; ++i)
    fmc->setTrans(0,i,i+1); // 0 is the initial state
  state=1;
  for(int i=0 ; i<=order ; ++i)
    {
      NthOrderStringIterator stateGen(i+1,alphabet);
      while(!stateGen.done())
	{
	  BOOM::String stateString=stateGen.getNextString();
	  BOOM::String history=
	    (stateString.length()<=order ? 
	     stateString : // at beginning of substrate; no drop necessary
	     stateString.substring(1,order));// drop one off to make room
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      BOOM::String nextStateString=history+alphabet.lookup(s);
	      int nextState=stringToState(nextStateString,*fmc);
	      fmc->setTrans(state,s,nextState);
	    }
	  ++state;
	}
    }  
  return fmc;
}



FastMarkovChain *MarkovChainCompiler::compileReverse(ContentSensor &mc)
{
  int order=mc.getOrder();
  FastMarkovChain *fmc=new FastMarkovChain(order,mc.getContentType());
  int alphabetSize=alphabet.getNumElements();

  // ===========================
  // Install state probabilities
  // ===========================
  Sequence seq;//dummy
  Symbol s;//dummy
  char c;//dummy
  fmc->setProb(0,LOG_OF_1);
  int state=1;
  for(int i=0 ; i<=order ; ++i)
    {
      NthOrderStringIterator stateGen(i+1,alphabet);
      while(!stateGen.done())
	{
	  BOOM::String stateString=stateGen.getNextString();
	  double score=mc.scoreSingleBase(seq,stateString,0,s,c);
	  fmc->setProb(state,score);
	  ++state;
	}
    }  

  // ===================
  // Install transitions
  // ===================

  for(int i=0 ; i<alphabetSize ; ++i)
    fmc->setTrans(0,i,0); // 0 is the final state, can never leave it
  state=1;
  for(int i=0 ; i<=order ; ++i)
    {
      NthOrderStringIterator stateGen(i+1,alphabet);
      while(!stateGen.done())
	{
	  BOOM::String stateString=stateGen.getNextString();
	  if(stateString.length()>order)
	    {
	      BOOM::String history=
		stateString.substring(1,order);// drop one off to make room
	      for(Symbol s=0 ; s<alphabetSize ; ++s)
		{
		  BOOM::String nextStateString=history+alphabet.lookup(s);
		  int nextState=stringToState(nextStateString,*fmc);
		  fmc->setTrans(state,s,nextState);
		}
	    }
	  ++state;
	}
    }  
  return fmc;
}



int MarkovChainCompiler::stringToState(const BOOM::String &s,
				       FastMarkovChain &fmc)
{
  int len=s.length();
  int state=1;
  for(int i=0 ; i<len-1 ; ++i)
    state+=fmc.statesOfOrder(i);
  int alphabetSize=alphabet.getNumElements();
  int factor=1;
  for(int i=0 ; i<len ; ++i)
    {
      state+=(alphabet.lookup(s[len-i-1]) * factor);
      factor*=alphabetSize;
    }
  return state;
}



BOOM::String MarkovChainCompiler::stateToString(int state,
					      FastMarkovChain &fmc)
{
  int alphabetSize=alphabet.getNumElements();
  int order=0;
  int firstStateOfOrder=1;
  for(; 1 ; ++order)
    {
      if(state<firstStateOfOrder) break;
      firstStateOfOrder+=fmc.statesOfOrder(order);
    }
  --order;
  firstStateOfOrder=1;
  for(int i=0 ; i<order ; ++i)
    firstStateOfOrder+=fmc.statesOfOrder(i);
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


