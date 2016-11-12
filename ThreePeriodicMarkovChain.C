/****************************************************************
 ThreePeriodicMarkovChain.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "ThreePeriodicMarkovChain.H"
#include <iostream>
#include <fstream>
#include "genezilla.H"
#include "Fast3PMC.H"


ThreePeriodicMarkovChain::ThreePeriodicMarkovChain(const BOOM::String &
						   filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file ")+filename
		   +" in ThreePeriodicMarkovChain()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="3P")
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      " into a ThreePeriodicMarkovChain (3P)";
  load(is);  
}



ThreePeriodicMarkovChain::ThreePeriodicMarkovChain(istream &is)
{
  load(is);
}



ThreePeriodicMarkovChain::ThreePeriodicMarkovChain(Strand strand,
						   ContentType contentType)
{
  setStrand(strand);
  setContentType(contentType);
}



ThreePeriodicMarkovChain::ThreePeriodicMarkovChain(BOOM::Vector<TrainingSequence*> &
						   sequences,int order,
						   int minSampleSize,
						   ContentType contentType)
{
  setContentType(contentType);
  setStrand(PLUS_STRAND);
  buildModels(sequences,minSampleSize,order);
}



ThreePeriodicMarkovChain::~ThreePeriodicMarkovChain()
{
  for(int i=0 ; i<3 ; ++i)
    delete chains[i];
}



double ThreePeriodicMarkovChain::scoreSingleBase(const Sequence &seq,
						 const BOOM::String &str,
						 int index,
						 Symbol s,char c)
{
  throw "noncoding version of scoreSingleBase() not applicable for 3PMC!";
}



void ThreePeriodicMarkovChain::scoreSingleBase(const Sequence &seq,
					       const BOOM::String &str,
					       int index,
					       Symbol s,char c,
					       double &scorePhase0,
					       double &scorePhase1,
					       double &scorePhase2)
{
  scorePhase0=chains[0]->scoreSingleBase(seq,str,index,s,c);
  scorePhase1=chains[1]->scoreSingleBase(seq,str,index,s,c);
  scorePhase2=chains[2]->scoreSingleBase(seq,str,index,s,c);
}



double ThreePeriodicMarkovChain::scoreSubsequence(const Sequence &seq,
						  const BOOM::String &str,
						  int begin,int length,
						  int seqPhase)
{
  // ### This is quick and dirty -- will be optimizing this class next week
  //     to use a suffix-tree-like data structure so that each base is
  //     read only once

  double score=0;
  int end=begin+length;

  switch(getStrand())
    {
    case FORWARD_STRAND:
      for(int pos=begin ; pos<end ; ++pos)
	score+=chains[(seqPhase+pos-begin)%3]->
	  scoreSingleBase(seq,str,pos,seq[pos],str[pos]);
      break;
    case REVERSE_STRAND:
      for(int pos=begin ; pos<end ; ++pos)
	score+=chains[posmod(seqPhase-(pos-begin))]->
	  scoreSingleBase(seq,str,pos,seq[pos],str[pos]);
      break;
    }
  return score;
}



bool ThreePeriodicMarkovChain::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good()) throw BOOM::String("Error creating file ")+filename+
		   "in ThreePeriodicMarkovChain::save()";
  return save(os);
}



void ThreePeriodicMarkovChain::buildModels(BOOM::Vector<TrainingSequence*> &sequences,
					   int minSampleSize,int order)
{
  for(int i=0 ; i<3 ; ++i)
    chains[i]=new MarkovChain(sequences,order,minSampleSize,i,
			      getContentType());
}



bool ThreePeriodicMarkovChain::save(ostream &os)
{
  os.precision(8);
  os << "3P" << endl;
  for(int i=0 ; i<3 ; ++i)
    chains[i]->save(os);
}



void ThreePeriodicMarkovChain::load(istream &is)
{
  BOOM::String dummy;
  for(int i=0 ; i<3 ; ++i)
    {
      is >> dummy;
      chains[i]=new MarkovChain(is);
    }
  ContentType contentType=chains[0]->getContentType();
  setContentType(contentType);
  setStrand(::getStrand(contentType));
}



ContentSensor *ThreePeriodicMarkovChain::reverseComplement()
{
  ThreePeriodicMarkovChain *other=
    new ThreePeriodicMarkovChain(complement(getStrand()),
				 ::reverseComplement(getContentType()));
  for(int i=0 ; i<3 ; ++i)
    other->chains[i]=
      static_cast<MarkovChain*>(chains[i]->reverseComplement());
  return other;
}



void ThreePeriodicMarkovChain::useLogOdds(ContentSensor &nullModel)
{
  ThreePeriodicMarkovChain &null3P=
    dynamic_cast<ThreePeriodicMarkovChain&>(nullModel);
  for(int i=0 ; i<3 ; ++i)
    chains[i]->useLogOdds(*null3P.chains[i]);
}



void ThreePeriodicMarkovChain::useLogOdds_anonymous(ContentSensor &nullModel)
{
  for(int i=0 ; i<3 ; ++i)
    chains[i]->useLogOdds_anonymous(nullModel);
}



ContentSensor *ThreePeriodicMarkovChain::compile()
{
  return new Fast3PMC(*this);
}


