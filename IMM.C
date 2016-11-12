/****************************************************************
 IMM.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "IMM.H"
#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/ProteinTrans.H"
#include "BOOM/Constants.H"
#include "BOOM/Alphabet.H"
#include "BOOM/Stacktrace.H"
#include "BOOM/Chi2FitTest.H"
#include "NthOrderStringIterator.H"
#include "MarkovChainCompiler.H"


extern Alphabet alphabet;


inline int hashTableSize(int order)
{
  // these are all primes

  switch(order)
    {
    case 0: return 11; break;
    case 1: return 29; break;
    case 2: return 127; break;
    case 3: return 619; break;
    default:
    case 4: return 3121; break;
    }
}



IMM::IMM(const IMM &other)
  : N(other.N), phase(other.phase), alphabetSize(other.alphabetSize),
    revComp(NULL), models(new BOOM::Vector<BOOM::StringMap<double>*>)
{
  for(int i=0 ; i<=N ; ++i)
    models->push_back(new BOOM::StringMap<double>(*(*other.models)[i]));
  setContentType(other.getContentType());
  setStrand(other.getStrand());
  if(getContentType()==INTERGENIC)
    revComp=this;
  else if(other.revComp && getStrand()==FORWARD_STRAND)
    revComp=new IMM(*other.revComp);
}



IMM::IMM(BOOM::Vector<TrainingSequence*> &v,int order,
	 int minSampleSize,int phase,ContentType contentType,
	 Strand strand)
  : N(order),
    alphabetSize(alphabet.getNumElements()),
    phase(phase),
    revComp(NULL),
    models(new BOOM::Vector<BOOM::StringMap<double>*>)
{
  setContentType(contentType);
  if(strand==EITHER_STRAND) strand=::getStrand(contentType);
  setStrand(strand);

  buildModels(v,minSampleSize);

  if(strand==FORWARD_STRAND)
    {
      BOOM::Vector<TrainingSequence*> rcSeqs;
      revCompSeqs(v,rcSeqs);
      revComp=new IMM(rcSeqs,order,minSampleSize,phase,
		      ::reverseComplement(contentType),
		      REVERSE_STRAND);
      revComp->revComp=this;
    }
}



IMM::IMM(istream &is,Strand strand)
  : revComp(NULL), models(new BOOM::Vector<BOOM::StringMap<double>*>)
{
  setStrand(strand);
  load(is);
}



IMM::IMM(const BOOM::String &filename)
  : revComp(NULL), models(new BOOM::Vector<BOOM::StringMap<double>*>)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file ")+filename
		   +" in IMM::IMM()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="IMM")
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      " into an IMM";
  load(is);
}



IMM::~IMM()
{
  int n=models->size();
  for(int i=0 ; i<n ; ++i)
    delete (*models)[i];
  delete models;

  /* ### THIS IS BAD!  THE REVCOMP WILL BE DELETED MULTIPLE TIMES...
         (UNLESS FORWARD_STRAND_ONLY IS DEFINED)
  if(getStrand()==FORWARD_STRAND && revComp!=this)
    delete revComp;
    */
}



double IMM::scoreSubsequence(const Sequence &seq,const BOOM::String &str,
				     int begin,int length,int)
{
  // This is slow -- please use FastIMM if you care about speed...

  Symbol dummySymbol;
  char dummyChar;
  double score=0;
  int end=begin+length;
  for(int pos=begin ; pos<end ; ++pos)
    score+=scoreSingleBase(seq,str,pos,dummySymbol,dummyChar);
  return score;
}



double IMM::scoreSingleBase(const Sequence &seq,const BOOM::String &str,
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
	    BOOM::StringMap<double> &model=*(*models)[order];
	    if(model.isDefined(p,index-order,order+1))
	      return model.lookup(p,index-order,order+1);
	  }
	throw BOOM::String("IMM::scoreSingleBase('+',")+
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
	    BOOM::StringMap<double> &model=*(*models)[order];
	    if(model.isDefined(p,index,order+1)) 
	      return model.lookup(p,index,order+1);
	  }
	throw BOOM::Stacktrace(
          BOOM::String("IMM::scoreSingleBase('-',")+
	    index+",strlen="+strlen(p)+",str="+
	  str.substring(index,maxOrder)+")");
      }

    default: throw BOOM::String(__FILE__)+__LINE__;
    }
}



void IMM::scoreSingleBase(const Sequence &seq,const BOOM::String &str,
				  int index,Symbol s,char c,
				  double &scorePhase0,double &scorePhase1,
				  double &scorePhase2)
{
  // I'm not a 3-periodic markov chain -- I'm just a plain, ordinary,
  // homogeneous markov chain:
  scorePhase0=scorePhase1=scorePhase2=scoreSingleBase(seq,str,index,s,c);
}



ContentSensor *IMM::reverseComplement()
{
  return revComp;
}



bool IMM::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good()) throw BOOM::String("Error creating file ")+filename+
		   "in IMM::save()";
  return save(os);
}



bool IMM::save(ostream &os)
{
  os.precision(8);
  os << "IMM" << endl;
  os << getContentType() << endl;
  os << N << "\t" << phase << "\t" << endl;
  int numModels=models->size();
  os << numModels << endl;
  for(int i=0 ; i<numModels ; ++i)
    {
      BOOM::StringMap<double> &model=*(*models)[i];
      os << model.size() << endl;
      BOOM::StringMap<double>::iterator cur=model.begin(), end=model.end();
      for(; cur!=end ; ++cur)
	os << (*cur).first << endl << (*cur).second << endl;
    }
  if(getStrand()==FORWARD_STRAND)  revComp->save(os);

  return true;
}



void IMM::load(istream &is)
{
  int numModels, numElements;
  BOOM::String str, pStr;
  ContentType contentType;
  
  is >> contentType >> N >> phase >> numModels;
  setContentType(contentType);
  setStrand(::getStrand(contentType)); // ### 9/17/2015

  for(int i=0 ; i<numModels ; ++i)
    {
      models->push_back(new BOOM::StringMap<double>(hashTableSize(N)));
      BOOM::StringMap<double> &model=*(*models)[i];
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
      revComp=new IMM(is,REVERSE_STRAND);
      revComp->revComp=this;
    }
}



void IMM::buildModels(BOOM::Vector<TrainingSequence*> &v,
			      int minSampleSize)
{
  /*
    This is the training procedure for IMMs
   */

  // Instantiate a separate hash table for each order
  counts=new BOOM::Vector<BOOM::StringMap<int>*>;
  for(int i=0 ; i<=N ; ++i)
    {
      int size=hashTableSize(N);
      models->push_back(new BOOM::StringMap<double>(size));
      counts->push_back(new BOOM::StringMap<int>(size));
    }

  // Install pseudocounts
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<int> &count=*(*counts)[order];
      NthOrderStringIterator iterator(order+1,alphabet);
      while(!iterator.done())
	{
	  BOOM::String ngram=iterator.getNextString();
	  count.lookup(ngram.c_str(),ngram.length())=1;
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
    {
      computeProbabilities_fw(minSampleSize);
      interpolate_fw();
    }
  else
    {
      computeProbabilities_rev(minSampleSize);
      interpolate_rev();
    }
  delete counts;
}



void IMM::computeProbabilities_fw(int minSampleSize)
{
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> &model=*(*models)[order];
      BOOM::StringMap<int> &count=*(*counts)[order];
      NthOrderStringIterator iterator(order,alphabet);
      while(!iterator.done())
	{
	  // Within this state (history), consider all emissions
	  BOOM::String history=iterator.getNextString();
	  int sampleSize=0;
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      BOOM::String ngram=history+alphabet.lookup(s);
	      if(count.isDefined(ngram.c_str(),ngram.length()))
		sampleSize+=int(count.lookup(ngram.c_str(),ngram.length()));
	    }
	  
	  // If sample size is insufficient, remove from the model
	  // (a lower-order model will always be available)
	  if(sampleSize<minSampleSize)
	    undefine_fw(history,count);

	  // Otherwise, normalize counts by sample size to produce a
	  // probability
	  else
	    for(Symbol s=0 ; s<alphabetSize ; ++s)
	      {
		BOOM::String ngram=history+alphabet.lookup(s);
		if(count.isDefined(ngram.c_str(),ngram.length()))
		  {
		    model.lookup(ngram.c_str(),ngram.length())=
		      log(count.lookup(ngram.c_str(),ngram.length())/
			  (double)sampleSize);
		  }
		else // ### this might be dangerous...
		  model.lookup(ngram.c_str(),ngram.length())=
		    NEGATIVE_INFINITY;
	      }
	}
    }
}



void IMM::computeProbabilities_rev(int minSampleSize)
{
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> &model=*(*models)[order];
      BOOM::StringMap<int> &count=*(*counts)[order];
      NthOrderStringIterator iterator(order,alphabet);
      while(!iterator.done())
	{
	  // Within this state (future), consider all emissions
	  BOOM::String future=iterator.getNextString();
	  int sampleSize=0;
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      BOOM::String ngram=alphabet.lookup(s)+future;
	      if(count.isDefined(ngram.c_str(),ngram.length()))
		sampleSize+=int(count.lookup(ngram.c_str(),ngram.length()));
	    }
	  
	  // If sample size is insufficient, remove from the model
	  // (a lower-order model will always be available)
	  if(sampleSize<minSampleSize)
	    undefine_rev(future,count);

	  // Otherwise, normalize counts by sample size to produce a
	  // probability
	  else
	    for(Symbol s=0 ; s<alphabetSize ; ++s)
	      {
		BOOM::String ngram=alphabet.lookup(s)+future;
		if(count.isDefined(ngram.c_str(),ngram.length()))
		  {
		    model.lookup(ngram.c_str(),ngram.length())=
		      log(count.lookup(ngram.c_str(),ngram.length())
			  /(double)sampleSize);
		  }
		else // ### this might be dangerous...
		  model.lookup(ngram.c_str(),ngram.length())=
		    NEGATIVE_INFINITY;
	      }
	}
    }
}



void IMM::undefine_fw(BOOM::String &history,
		      BOOM::StringMap<int> &model)
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



void IMM::undefine_rev(BOOM::String &future,
		       BOOM::StringMap<int> &model)
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



void IMM::updateCounts_fw(BOOM::String &str,int order,int seqPhase,
				  int boostCount)
{
  /*
    This method simply counts all the n-grams of length order+1.
   */

  BOOM::StringMap<int> &count=*(*counts)[order];
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
	throw BOOM::String("seq phase=-1 in IMM::updateCounts(), "
	  "called as part of ThreePeriodicIMM:\n")+str;
      firstPos=(phase-seqPhase+3)%3;
      periodicity=3;
    }

  for(int pos=firstPos ; pos<len ; pos+=periodicity)
    {
      int begin=pos-order;
      if(begin<0) continue;
      if(count.isDefined(str.c_str(),begin,windowSize)) 
	count.lookup(str.c_str(),begin,windowSize)+=boostCount;
      else 
	count.lookup(str.c_str(),begin,windowSize)=boostCount;
    }
}



void IMM::updateCounts_rev(BOOM::String &str,int order,int seqPhase,
				   int boostCount)
{
  /*
    This method simply counts all the n-grams of length order+1.
   */

  BOOM::StringMap<int> &count=*(*counts)[order];
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
	throw BOOM::String("seq phase=-1 in IMM::updateCounts(), "
	  "called as part of ThreePeriodicIMM (-):\n")+str;
      firstPos=(seqPhase-phase+3)%3;
      periodicity=3;
    }

  for(int pos=firstPos ; pos+order<len ; pos+=periodicity)
    {
      if(count.isDefined(str.c_str(),pos,windowSize)) 
	count.lookup(str.c_str(),pos,windowSize)+=boostCount;
      else 
	count.lookup(str.c_str(),pos,windowSize)=boostCount;
    }
}



void IMM::revCompSeqs(BOOM::Vector<TrainingSequence*> &forwardSeqs,
		      BOOM::Vector<TrainingSequence*> &revSeqs)
{
  BOOM::Vector<TrainingSequence*>::iterator cur=forwardSeqs.begin(), 
    end=forwardSeqs.end();
  for(; cur!=end ; ++cur)
    revSeqs.push_back(
      static_cast<TrainingSequence*>((*cur)->reverseComplement(alphabet)));
}



void IMM::useLogOdds_anonymous(ContentSensor &nullModel)
{
  Strand strand=getStrand();
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> *model=(*models)[order];
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
	  cout << "="<<nullScore<<" vs "<<(*cur).second <<p<<endl;
	  (*cur).second-=nullScore;
	}
    }

  if(strand==FORWARD_STRAND && revComp && revComp!=this)
    revComp->useLogOdds_anonymous(*nullModel.reverseComplement());
}



void IMM::useLogOdds(ContentSensor &nullModel)
{
  IMM &other=dynamic_cast<IMM&>(nullModel);
  for(int order=0 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> *model=(*models)[order];
      BOOM::StringMap<double> *otherModel=(*other.models)[order];
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



void IMM::interpolate_fw()
{
  BOOM::Vector<BOOM::StringMap<double>*> interpolated;
  interpolated.push_back((*models)[0]);
  
  for(int order=1 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> *newModel=
	new BOOM::StringMap<double>(hashTableSize(order));
      interpolated.push_back(newModel);
      BOOM::StringMap<double> *prevIMM=interpolated[order-1];
      BOOM::StringMap<double> *model=(*models)[order];
      BOOM::StringMap<int> *count=(*counts)[order];
      BOOM::StringMap<int> *smallerCount=(*counts)[order-1];
      NthOrderStringIterator iterator(order,alphabet);
      while(!iterator.done())
	{
	  BOOM::String longHistory=iterator.getNextString();
	  BOOM::String shortHistory=longHistory.substr(1,order-1);
	  BOOM::Vector<int> longCounts, shortCounts; // distributions
	  int longTotal=0;
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      char c=alphabet.lookup(s);
	      BOOM::String ngram=longHistory+c;
	      const char *p=ngram.c_str(); int len=ngram.length();
	      if(!count->isDefined(p,len)) continue;
	      int longCount=count->lookup(p,len);
	      longCounts.push_back(longCount);
	      longTotal+=longCount;
	      ngram=shortHistory+c;
	      int shortCount=
		smallerCount->lookup(ngram.c_str(),ngram.length());
	      /*
	      // ### the following is true to the original paper...
	      int shortCount=
		10000*prevIMM->lookup(ngram.c_str(),ngram.length());
	       */
	      shortCounts.push_back(shortCount);
	    }

	  double p=chiTest(longCounts,shortCounts);

	  double confidence=1-p;// confidence that the distributions differ
	  if(confidence<0.5) confidence=0.0;
	  double samplePortion=longTotal/double(400);
	  double lambda=
	    (samplePortion>=1 ? 1.0 : confidence*samplePortion);
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      char c=alphabet.lookup(s);
	      BOOM::String ngram=longHistory+c;
	      const char *p=ngram.c_str(); int len=ngram.length();
	      if(!count->isDefined(p,len)) continue;
	      newModel->lookup(p,len)=
		lambda * model->lookup(p,len) +
		(1-lambda) * prevIMM->lookup(p+1,len-1);
	    }
	}
    }
  for(int i=1 ; i<=N ; ++i) delete (*models)[i];
  *models=interpolated;
}



double IMM::chiTest(BOOM::Vector<int> &counts1,BOOM::Vector<int> &counts2)
{
  int sum1=0, sum2=0, n1=counts1.size(), n2=counts2.size();
  for(int i=0 ; i<n1 ; ++i) sum1+=counts1[i];
  for(int i=0 ; i<n2 ; ++i) sum2+=counts2[i];
  BOOM::Vector<int> expectedCounts;
  double denom=sum2;
  for(int i=0 ; i<n2 ; ++i) 
    expectedCounts.push_back(int(counts2[i]/denom*sum1+0.5));
  BOOM::Chi2FitTest test(counts1,expectedCounts,chiSquaredTable);
  return test.getP();
}



void IMM::interpolate_rev()
{
  BOOM::Vector<BOOM::StringMap<double>*> interpolated;
  interpolated.push_back((*models)[0]);
  
  for(int order=1 ; order<=N ; ++order)
    {
      BOOM::StringMap<double> *newModel=
	new BOOM::StringMap<double>(hashTableSize(order));
      interpolated.push_back(newModel);
      BOOM::StringMap<double> *prevIMM=interpolated[order-1];
      BOOM::StringMap<double> *model=(*models)[order];
      BOOM::StringMap<int> *count=(*counts)[order];
      BOOM::StringMap<int> *smallerCount=(*counts)[order-1];
      NthOrderStringIterator iterator(order,alphabet);
      while(!iterator.done())
	{
	  BOOM::String longFuture=iterator.getNextString();
	  BOOM::String shortFuture=longFuture.substr(0,order-1);
	  BOOM::Vector<int> longCounts, shortCounts; // distributions
	  int longTotal=0;
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      char c=alphabet.lookup(s);
	      BOOM::String ngram=c+longFuture;
	      const char *p=ngram.c_str(); int len=ngram.length();
	      if(!count->isDefined(p,len)) continue;
	      int longCount=count->lookup(p,len);
	      longCounts.push_back(longCount);
	      longTotal+=longCount;
	      ngram=c+shortFuture;
	      int shortCount=
		smallerCount->lookup(ngram.c_str(),ngram.length());
	      shortCounts.push_back(shortCount);
	    }
	  
	  double p=chiTest(longCounts,shortCounts);

	  double confidence=1-p;// confidence that the distributions differ
	  if(confidence<0.5) confidence=0.0;
	  double samplePortion=longTotal/double(400);
	  double lambda=
	    (samplePortion>=1 ? 1.0 : confidence*samplePortion);
	  for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      char c=alphabet.lookup(s);
	      BOOM::String ngram=c+longFuture;
	      const char *p=ngram.c_str(); int len=ngram.length();
	      if(!count->isDefined(p,len)) continue;
	      newModel->lookup(p,len)=
		lambda * model->lookup(p,len) +
		(1-lambda) * prevIMM->lookup(p,len-1);
	    }
	}
    }

  for(int i=1 ; i<=N ; ++i) delete (*models)[i];
  *models=interpolated;
}



ContentSensor *IMM::compile()
{
  return MarkovChainCompiler::compile(*this);
}

