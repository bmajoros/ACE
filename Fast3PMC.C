/****************************************************************
 Fast3PMC.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "Fast3PMC.H"
#include <iostream>
#include "MarkovChainCompiler.H"


Fast3PMC::Fast3PMC(ThreePeriodicMarkovChain &slowModel)
  : revComp(NULL)
{
  ContentType contentType=slowModel.getContentType();
  setContentType(contentType);
  setStrand(::getStrand(contentType));
  compileFrom(slowModel);
}



Fast3PMC::Fast3PMC(ThreePeriodicIMM &slowModel)
  : revComp(NULL)
{
  ContentType contentType=slowModel.getContentType();
  setContentType(contentType);
  setStrand(::getStrand(contentType));
  compileFrom(slowModel);
}



Fast3PMC::Fast3PMC(ContentType contentType)
  : revComp(NULL)
{
  setContentType(contentType);
  setStrand(::getStrand(contentType));
  chains[0]=chains[1]=chains[2]=NULL;
}



Fast3PMC::Fast3PMC(const BOOM::String &filename)
  : revComp(NULL)
{
  load(filename);
}



Fast3PMC::Fast3PMC(BOOM::File &file)
  : revComp(NULL)
{
  load(file);
}



double Fast3PMC::scoreSingleBase(const Sequence &seq,
				 const BOOM::String &str,
				 int index,Symbol s,char c)
{
  throw "attempt to use Fast3PMC on a non-periodic feature";
}



void Fast3PMC::scoreSingleBase(const Sequence &seq,
			       const BOOM::String &str,int index,
			       Symbol s,char c,double &scorePhase0,
			       double &scorePhase1,double &scorePhase2)
{
  scorePhase0=chains[0]->scoreSingleBase(seq,str,index,s,c);
  scorePhase1=chains[1]->scoreSingleBase(seq,str,index,s,c);
  scorePhase2=chains[2]->scoreSingleBase(seq,str,index,s,c);
}



double Fast3PMC::scoreSubsequence(const Sequence &seq,
				  const BOOM::String &str,
				  int begin,int length,int phase)
{
  char c;
  double score=0;
  int end=begin+length;
  double s0, s1, s2;
  for(int pos=begin ; pos<end ; ++pos)
    {
      scoreSingleBase(seq,str,pos,seq[pos],c,s0,s1,s2);
      switch(phase)
	{
	case 0: score+=s0; break;
	case 1: score+=s1; break;
	case 2: score+=s2; break;
	}
      phase=(phase+1)%3;
    }
  return score;
}



ContentSensor *Fast3PMC::reverseComplement()
{
  if(!revComp)
    {
      revComp=new Fast3PMC(::reverseComplement(getContentType()));
      revComp->revComp=this;
      for(int i=0 ; i<3 ; ++i)
	revComp->chains[i]=
	  static_cast<FastMarkovChain*>(chains[i]->reverseComplement());
    }
  return revComp;
}



bool Fast3PMC::save(const BOOM::String &filename)
{
  BOOM::File file(filename,"w");
  return save(file);
}



bool Fast3PMC::save(ostream &os)
{
  throw "Attempt to save Fast3PMC into ASCII file -- please use binary files";
}



bool Fast3PMC::save(BOOM::File &file)
{
  file << "Fast3PMC" << static_cast<int>(getContentType());
  for(int i=0 ; i<3 ; ++i)
    if(!chains[i]->save(file)) return false;
  return true;
}



void Fast3PMC::load(const BOOM::String &filename)
{
  BOOM::File file(filename);
  load(file);
}



void Fast3PMC::load(BOOM::File &file)
{
  BOOM::String modelType;
  file >> modelType;
  if(modelType!="Fast3PMC")
    throw BOOM::String("Attempt to load a model of type ")+modelType+
      " into Fast3PMC (optimized 3-period markov chain)";
  int contentType;
  file >> contentType;
  setContentType(static_cast<ContentType>(contentType));
  setStrand(::getStrand(static_cast<ContentType>(contentType)));
  for(int i=0 ; i<3 ; ++i)
    chains[i]=new FastMarkovChain(file);
}



void Fast3PMC::reset(const Sequence &seq,const BOOM::String &str,int pos)
{
  for(int i=0 ; i<3 ; ++i)
    chains[i]->reset(seq,str,pos);
}



void Fast3PMC::compileFrom(ThreePeriodicMarkovChain &slowModel)
{
  int order=slowModel.getOrder();
  for(int i=0 ; i<3 ; ++i)
    chains[i]=MarkovChainCompiler::compile(*slowModel.chains[i]);
}



void Fast3PMC::compileFrom(ThreePeriodicIMM &slowModel)
{
  int order=slowModel.getOrder();
  for(int i=0 ; i<3 ; ++i)
    chains[i]=MarkovChainCompiler::compile(*slowModel.chains[i]);
}


