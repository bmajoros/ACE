/****************************************************************
 ThreePeriodicIMM.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "ThreePeriodicIMM.H"
#include <iostream>
#include <fstream>
#include "genezilla.H"
#include "Fast3PMC.H"


ThreePeriodicIMM::ThreePeriodicIMM(const BOOM::String &filename)
  : revComp(NULL)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file ")+filename
		   +" in ThreePeriodicIMM()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="3PIMM")
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      " into a ThreePeriodicIMM (3PIMM)";
  load(is);  
}



ThreePeriodicIMM::ThreePeriodicIMM(istream &is)
  : revComp(NULL)
{
  load(is);
}



ThreePeriodicIMM::ThreePeriodicIMM(Strand strand,ContentType contentType)
  : revComp(NULL)
{
  setStrand(strand);
  setContentType(contentType);
}



ThreePeriodicIMM::ThreePeriodicIMM(BOOM::Vector<TrainingSequence*> &
				   sequences,int order,
				   int minSampleSize,
				   ContentType contentType)
  : revComp(NULL)
{
  setContentType(contentType);
  setStrand(PLUS_STRAND);
  buildModels(sequences,minSampleSize,order);
}



ThreePeriodicIMM::~ThreePeriodicIMM()
{
  for(int i=0 ; i<3 ; ++i)
    delete chains[i];
}



double ThreePeriodicIMM::scoreSingleBase(const Sequence &seq,
						 const BOOM::String &str,
						 int index,
						 Symbol s,char c)
{
  throw "noncoding version of scoreSingleBase() not applicable for 3PIMM!";
}



void ThreePeriodicIMM::scoreSingleBase(const Sequence &seq,
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



double ThreePeriodicIMM::scoreSubsequence(const Sequence &seq,
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



bool ThreePeriodicIMM::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good()) throw BOOM::String("Error creating file ")+filename+
		   "in ThreePeriodicIMM::save()";
  return save(os);
}



void ThreePeriodicIMM::buildModels(BOOM::Vector<TrainingSequence*> &sequences,
					   int minSampleSize,int order)
{
  for(int i=0 ; i<3 ; ++i)
    chains[i]=new IMM(sequences,order,minSampleSize,i,
			      getContentType());
}



bool ThreePeriodicIMM::save(ostream &os)
{
  os.precision(8);
  os << "3PIMM" << endl;
  for(int i=0 ; i<3 ; ++i)
    chains[i]->save(os);
}



void ThreePeriodicIMM::load(istream &is)
{
  BOOM::String dummy;
  for(int i=0 ; i<3 ; ++i)
    {
      is >> dummy;
      chains[i]=new IMM(is);
    }
  ContentType contentType=chains[0]->getContentType();
  setContentType(contentType);
  setStrand(::getStrand(contentType));
}



ContentSensor *ThreePeriodicIMM::reverseComplement()
{
  if(!revComp) {
    revComp=new ThreePeriodicIMM(complement(getStrand()),
				 ::reverseComplement(getContentType()));
    for(int i=0 ; i<3 ; ++i)
      revComp->chains[i]=
	static_cast<IMM*>(chains[i]->reverseComplement());
  }
  return revComp;
}



void ThreePeriodicIMM::useLogOdds(ContentSensor &nullModel)
{
  ThreePeriodicIMM &null3P=
    dynamic_cast<ThreePeriodicIMM&>(nullModel);
  for(int i=0 ; i<3 ; ++i)
    chains[i]->useLogOdds(*null3P.chains[i]);
}



void ThreePeriodicIMM::useLogOdds_anonymous(ContentSensor &nullModel)
{
  for(int i=0 ; i<3 ; ++i)
    chains[i]->useLogOdds_anonymous(nullModel);
}



ContentSensor *ThreePeriodicIMM::compile()
{
  return new Fast3PMC(*this);
}
