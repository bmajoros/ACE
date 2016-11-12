/****************************************************************
 MarkovChain.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "MarkovChain.H"
#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/ProteinTrans.H"
#include "BOOM/Constants.H"
#include "BOOM/Alphabet.H"
#include "BOOM/Stacktrace.H"
#include "NthOrderStringIterator.H"
#include "MarkovChainCompiler.H"


extern Alphabet alphabet;


inline int hashTableSize(int order)
{
  // these are all primes

  switch(order)
    {
    case 0: return 11;
    case 1: return 53/*29*/;
    case 2: return 293/*127*/;
    case 3: return 1201/*619*/;
    case 4: return 3121;//6199/*3121*/;
    case 5: return 3121;//30983;
    default:
    case 6: return 3121;//199999;
    }
}



MarkovChain::MarkovChain(const MarkovChain &other)
  : N(other.N), phase(other.phase), alphabetSize(other.alphabetSize),
    revComp(NULL)
{
  for(int i=0 ; i<=N ; ++i)
    models.push_back(new BOOM::StringMap<double>(*other.models[i]));
  setContentType(other.getContentType());
  setStrand(other.getStrand());
  if(getContentType()==INTERGENIC)
    revComp=this;
  else if(other.revComp && getStrand()==FORWARD_STRAND)
    revComp=new MarkovChain(*other.revComp);
}



MarkovChain::MarkovChain(BOOM::Vector<TrainingSequence*> &v,int order,
			 int minSampleSize,int phase,
			 ContentType contentType,
			 Strand strand)
  : N(order),
    alphabetSize(alphabet.getNumElements()),
    phase(phase),
    revComp(NULL)
{
  setContentType(contentType);
  if(strand==EITHER_STRAND) strand=::getStrand(contentType);
  setStrand(strand);

  buildModels(v,minSampleSize);

  if(strand==FORWARD_STRAND)
    {
      BOOM::Vector<TrainingSequence*> rcSeqs;
      revCompSeqs(v,rcSeqs);
      revComp=new MarkovChain(rcSeqs,order,minSampleSize,phase,
			      ::reverseComplement(contentType),
			      REVERSE_STRAND);
      revComp->revComp=this;
    }
}



MarkovChain::MarkovChain(istream &is,Strand strand)
  : revComp(NULL)
{
  setStrand(strand);
  load(is);
}



MarkovChain::MarkovChain(const BOOM::String &filename)
  : revComp(NULL)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file ")+filename
		   +" in MC::MC()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="MC")
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      " into an MC (MarkovChain)";
  setStrand(FORWARD_STRAND);
  load(is);
}



MarkovChain::~MarkovChain()
{
  int n=models.size();
  for(int i=0 ; i<n ; ++i)
    delete models[i];
}



double MarkovChain::scoreSubsequence(const Sequence &seq,
				     const BOOM::String &str,
				     int begin,int length,int)
{
  // This is slow -- please use FastMarkovChain if you care about speed...

  Symbol dummySymbol;
  char dummyChar;
  double score=0;
  int end=begin+length;
  for(int pos=begin ; pos<end ; ++pos)
    score+=scoreSingleBase(seq,str,pos,dummySymbol,dummyChar);
  return score;
}



double MarkovChain::scoreSingleBase(const Sequence &seq,
				    const BOOM::String &str,
				    int index,Symbol s,char c)
{
  const char *p=str.c_str();
  switch(getStrand())
    {
    case PLUS_STRAND:
      {
	int maxOrder=(index>N ? N : index);
	for(int order=maxOrder ; order>=0 ; --order)
	  {
	    BOOM::StringMap<double> &model=*models[order];
	    if(model.isDefined(p,index-order,order+1))
	      {
		double d=model.lookup(p,index-order,order+1);
		//cout<<index<<"="<<d<<endl;//###
		return d;
	      }
	  }
	throw BOOM::String("MarkovChain::scoreSingleBase('+',")+
	  index+",strlen="+strlen(p)+",str="+
	  str.substring(index,maxOrder)+")";
      }

    case MINUS_STRAND:
      {
	/*
	  On the minus strand we have to take our contexts from the
	  right (but only because we trained the model that way)
	 */
	int seqLen=str.length();
	int maxOrder=seqLen-index-1;
	if(maxOrder>N) maxOrder=N;
	for(int order=maxOrder ; order>=0 ; --order)
	  {
	    BOOM::StringMap<double> &model=*models[order];
	    if(model.isDefined(p,index,order+1)) 
	      {
		double d=model.lookup(p,index,order+1);
		//cout<<index<<"="<<d<<" "<<getStrand()<<endl;//###
		return d;
	      }
	  }
	throw BOOM::Stacktrace(
          BOOM::String("MarkovChain::scoreSingleBase('-',")+
	    index+",strlen="+strlen(p)+",str="+
	  str.substring(index,maxOrder)+")");
      }

    default: throw BOOM::String(__FILE__)+" : line "+__LINE__;
    }
}



void MarkovChain::scoreSingleBase(const Sequence &seq,const BOOM::String &str,
				  int index,Symbol s,char c,
				  double &scorePhase0,double &scorePhase1,
				  double &scorePhase2)
{
  // I'm not a 3-periodic markov chain -- I'm just a plain, ordinary,
  // homogeneous markov chain:
  scorePhase0=scorePhase1=scorePhase2=scoreSingleBase(seq,str,index,s,c);
}



ContentSensor *MarkovChain::reverseComplement()
{
  return revComp;
}



bool MarkovChain::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good()) throw BOOM::String("Error creating file ")+filename+
		   "in MarkovChain::save()";
  return save(os);
}



bool MarkovChain::save(ostream &os)
{
  os.precision(8);
  os << "MC" << endl;
  os << getContentType() << endl;
  os << N << "\t" << phase << "\t" << endl;
  int numModels=models.size();
  os << numModels << endl;
  for(int i=0 ; i<numModels ; ++i)
    {
      BOOM::StringMap<double> &model=*models[i];
      os << model.size() << endl;
      BOOM::StringMap<double>::iterator cur=model.begin(), end=model.end();
      for(; cur!=end ; ++cur)
	os << (*cur).first << endl << (*cur).second << endl;
    }
  if(getStrand()==FORWARD_STRAND)
    revComp->save(os);
  return true;
}



void MarkovChain::load(istream &is)
{
  int numModels, numElements;
  BOOM::String str, pStr;
  ContentType contentType;
  
  is >> contentType >> N >> phase >> numModels;
  setContentType(contentType);

  for(int i=0 ; i<numModels ; ++i)
    {
      models.push_back(new BOOM::StringMap<double>(hashTableSize(N)));
      BOOM::StringMap<double> &model=*models[i];
      is >> numElements;
      for(int j=0 ; j<numElements ; ++j)
	{
	  is >> str >> pStr;
	  model.lookup(str.c_str(),str.length())=pStr.asDouble();
	}
    }

  if(getStrand()==FORWARD_STRAND)
    {
      BOOM::String modelType;
      is >> modelType;
      revComp=new MarkovChain(is,REVERSE_STRAND);
      revComp->revComp=this;
    }
}



void MarkovChain::buildModels(BOOM::Vector<TrainingSequence*> &v,
			      int minSampleSize)
{
  /*
    This is the training procedure for Markov chains.
   */

  // Instantiate a separate hash table for each order
  for(int i=0 ; i<=N ; ++i)
    models.push_back(new BOOM::StringMap<double>(hashTableSize(N)));

  // Install pseudocounts
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> &model=*models[order];
      NthOrderStringIterator iterator(order+1,alphabet);
      while(!iterator.done())
	{
	  BOOM::String ngram=iterator.getNextString();
	  model.lookup(ngram.c_str(),ngram.length())=1;
	}
    }

  // Extract all windows up to size N and count the different
  // combinations that actually occur
  BOOM::Vector<TrainingSequence*>::iterator cur=v.begin(), end=v.end();
  for(; cur!=end ; ++cur)
    {
      TrainingSequence &seq=**cur;
      BOOM::String &str=*seq.toString(alphabet);
      for(int order=0 ; order<=N ; ++order)
	if(getStrand()==FORWARD_STRAND) 
	  updateCounts_fw(str,order,seq.getPhase(),seq.getBoostCount());
	else
	  updateCounts_rev(str,order,seq.getPhase(),seq.getBoostCount());
      delete &str;
    }
  
  // Now convert the counts into conditional probabilities
  if(getStrand()==FORWARD_STRAND) 
    computeProbabilities_fw(minSampleSize);
  else 
    computeProbabilities_rev(minSampleSize);
}



void MarkovChain::computeProbabilities_fw(int minSampleSize)
{
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> &model=*models[order];
      NthOrderStringIterator iterator(order,alphabet);
      while(!iterator.done())
	{
	  // Within this state (history), consider all emissions
	  BOOM::String history=iterator.getNextString();
	  int sampleSize=0;
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      BOOM::String ngram=history+alphabet.lookup(s);
	      if(model.isDefined(ngram.c_str(),ngram.length()))
		sampleSize+=int(model.lookup(ngram.c_str(),ngram.length()));
	    }
	  
	  // If sample size is insufficient, remove from the model
	  // (a lower-order model will always be available)
	  if(sampleSize<minSampleSize)
	    undefine_fw(history,model);

	  // Otherwise, normalize counts by sample size to produce a
	  // probability
	  else
	    for(Symbol s=0 ; s<alphabetSize ; ++s)
	      {
		BOOM::String ngram=history+alphabet.lookup(s);
		if(model.isDefined(ngram.c_str(),ngram.length()))
		  {
		    double &d=model.lookup(ngram.c_str(),ngram.length());
		    d=log(d/sampleSize);
		  }
		else // ### this might be dangerous...
		  model.lookup(ngram.c_str(),ngram.length())=
		    NEGATIVE_INFINITY;
	      }
	}
    }
}



void MarkovChain::computeProbabilities_rev(int minSampleSize)
{
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> &model=*models[order];
      NthOrderStringIterator iterator(order,alphabet);
      while(!iterator.done())
	{
	  // Within this state (future), consider all emissions
	  BOOM::String future=iterator.getNextString();
	  int sampleSize=0;
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      BOOM::String ngram=alphabet.lookup(s)+future;
	      if(model.isDefined(ngram.c_str(),ngram.length()))
		sampleSize+=int(model.lookup(ngram.c_str(),ngram.length()));
	    }
	  
	  // If sample size is insufficient, remove from the model
	  // (a lower-order model will always be available)
	  if(sampleSize<minSampleSize)
	    undefine_rev(future,model);

	  // Otherwise, normalize counts by sample size to produce a
	  // probability
	  else
	    for(Symbol s=0 ; s<alphabetSize ; ++s)
	      {
		BOOM::String ngram=alphabet.lookup(s)+future;
		if(model.isDefined(ngram.c_str(),ngram.length()))
		  {
		    double &d=model.lookup(ngram.c_str(),ngram.length());
		    d=log(d/sampleSize);
		  }
		else // ### this might be dangerous...
		  model.lookup(ngram.c_str(),ngram.length())=
		    NEGATIVE_INFINITY;
	      }
	}
    }
}



void MarkovChain::undefine_fw(BOOM::String &history,
			      BOOM::StringMap<double> &model)
{
  /*
    This method removes from the given model all n-grams consisting
    of the given history followed by a single base.
   */

  for(Symbol s=0 ; s<alphabetSize ; ++s)
    {
      BOOM::String ngram=history+alphabet.lookup(s);
      if(model.isDefined(ngram.c_str(),ngram.length())) 
	model.remove(ngram.c_str(),ngram.length());
    }
}



void MarkovChain::undefine_rev(BOOM::String &future,
			      BOOM::StringMap<double> &model)
{
  /*
    This method removes from the given model all n-grams consisting
    of the given future preceded by a single base.
   */

  for(Symbol s=0 ; s<alphabetSize ; ++s)
    {
      BOOM::String ngram=alphabet.lookup(s)+future;
      if(model.isDefined(ngram.c_str(),ngram.length())) 
	model.remove(ngram.c_str(),ngram.length());
    }
}



void MarkovChain::updateCounts_fw(BOOM::String &str,int order,int seqPhase,
				  int boostCount)
{
  /*
    This method simply counts all the n-grams of length order+1.
   */

  BOOM::StringMap<double> &model=*models[order];
  int windowSize=order+1;
  int len=str.length();

  int firstPos, periodicity;
  if(phase==NO_PHASE) 
    {
      firstPos=0;
      periodicity=1;
    }
  else 
    {
      if(seqPhase==-1) 
	throw BOOM::String("seq phase=-1 in MarkovChain::updateCounts(), "
	  "called as part of ThreePeriodicMarkovChain:\n")+str;
      firstPos=(phase-seqPhase+3)%3;
      periodicity=3;
    }

  for(int pos=firstPos ; pos<len ; pos+=periodicity)
    {
      int begin=pos-order;
      if(begin<0) continue;
      if(model.isDefined(str.c_str(),begin,windowSize)) 
	model.lookup(str.c_str(),begin,windowSize)+=boostCount;
      else 
	model.lookup(str.c_str(),begin,windowSize)=boostCount;
    }
}



void MarkovChain::updateCounts_rev(BOOM::String &str,int order,int seqPhase,
				   int boostCount)
{
  /*
    This method simply counts all the n-grams of length order+1.
   */

  BOOM::StringMap<double> &model=*models[order];
  int windowSize=order+1;
  int len=str.length();

  int firstPos, periodicity;
  if(phase==NO_PHASE) 
    {
      firstPos=0;
      periodicity=1;
    }
  else 
    {
      if(seqPhase==-1) 
	throw BOOM::String("seq phase=-1 in MarkovChain::updateCounts(), "
	  "called as part of ThreePeriodicMarkovChain (-):\n")+str;
      firstPos=(seqPhase-phase+3)%3;
      periodicity=3;
    }

  for(int pos=firstPos ; pos+order<len ; pos+=periodicity)
    {
      if(model.isDefined(str.c_str(),pos,windowSize)) 
	model.lookup(str.c_str(),pos,windowSize)+=boostCount;
      else 
	model.lookup(str.c_str(),pos,windowSize)=boostCount;
    }
}



void MarkovChain::revCompSeqs(BOOM::Vector<TrainingSequence*> &forwardSeqs,
			      BOOM::Vector<TrainingSequence*> &revSeqs)
{
  BOOM::Vector<TrainingSequence*>::iterator cur=forwardSeqs.begin(), 
    end=forwardSeqs.end();
  for(; cur!=end ; ++cur)
    revSeqs.push_back(
      static_cast<TrainingSequence*>((*cur)->reverseComplement(alphabet)));
}



void MarkovChain::useLogOdds_anonymous(ContentSensor &nullModel)
{
  Strand strand=getStrand();
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> *model=models[order];
      BOOM::StringMap<double>::iterator cur=model->begin(),
	end=model->end();
      for(; cur!=end ; ++cur)
	{
	  const char *p=(*cur).first;
	  int len=strlen(p);
	  Sequence seq(p,alphabet);
	  int targetPos=(strand==FORWARD_STRAND ? len-1 : 0);
	  double nullScore=nullModel.scoreSingleBase(seq,p,targetPos,
						     seq[targetPos],
						     p[targetPos]);
	  (*cur).second-=nullScore;
	}
    }

  if(strand==FORWARD_STRAND && revComp && revComp!=this)
    revComp->useLogOdds_anonymous(*nullModel.reverseComplement());
}



void MarkovChain::useLogOdds(ContentSensor &nullModel)
{
  MarkovChain &other=dynamic_cast<MarkovChain&>(nullModel);
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> *model=models[order];
      BOOM::StringMap<double> *otherModel=other.models[order];
      BOOM::StringMap<double>::iterator cur=model->begin(),
	end=model->end();
      for(; cur!=end ; ++cur)
	{
	  const char *p=(*cur).first;
	  int len=strlen(p);
	  (*cur).second-=otherModel->lookup(p,0,len);
	}
    }

  if(getStrand()==FORWARD_STRAND && revComp && revComp!=this)
    revComp->useLogOdds(*nullModel.reverseComplement());
}



ContentSensor *MarkovChain::compile()
{
  return MarkovChainCompiler::compile(*this);
}

