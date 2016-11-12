/**************************************************************
 WAM.C
 Copyright (C)2003 Mihaela Pertea (mpertea@tigr.org)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include "WAM.H"
#include <iostream>
#include <fstream>


WAM::WAM(GarbageCollector &gc,const WAM &other,bool revComp)  
  : SignalSensor(gc,other,revComp) 
{
  // copy ctor
  
  if(!revComp) {
    int n=other.matrix.size();  
    for(int i=0;i<n; ++i)
      matrix.push_back(new MarkovChain(*other.matrix[i]));

  }
}



WAM::WAM(GarbageCollector &gc,BOOM::String &filename)
  : SignalSensor(gc)
{
  // ctor

  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file ")+filename+
		   "in WAM::WAM()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="WAM") 
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      "in into a WAM";
  load(is);
}



WAM::WAM(GarbageCollector &gc,istream &is)  
  : SignalSensor(gc)
{
  // ctor

  load(is);
}



WAM::WAM(GarbageCollector &gc,BOOM::Vector<TrainingSequence*> &sequences,
	 int order,int minSampleSize,SignalType signalType,
	 int consensusOffset,int consensusLength)
  : SignalSensor(gc)
{
  // ctor

  BOOM::Vector<TrainingSequence*> examples;

  setStrand(FORWARD_STRAND);
  setSignalType(signalType);
  int windowSize=sequences[0]->getLength() ;
  setSizes(consensusLength,consensusOffset,windowSize);   

  int nosequences=sequences.size();

  for(int i=0; i<windowSize; ++i) 
    {
      // just extract the necessary examples from the &sequences
      
      int min = 0 > i-order ? 0 : i-order;
      
      // ### this considers 0-order markov chains for the first base 
      // in the window; we might want to try higher orders for this
      
      for(int j=0;j<nosequences;j++) 
	{
	  TrainingSequence *s=new TrainingSequence();
	  sequences[j]->getSubsequence(min,i-min+1,*s);
	  examples.push_back(s);
	}
      
      matrix.push_back(new MarkovChain(examples,order,minSampleSize,-1,
				       UNKNOWN_CONTENT_FORWARD));
      
      BOOM::Vector<TrainingSequence*>::iterator cur=examples.begin(), 
	end=examples.end();
      for(; cur!=end ; ++cur) delete *cur;
      examples.clear();
    }
}



WAM::~WAM()
{
  BOOM::Vector<MarkovChain *>::iterator cur=matrix.begin(), end=matrix.end();
  for(; cur!=end ; ++cur)
    delete *cur;
}



double WAM::getLogP(const Sequence &seq,const BOOM::String &str,int begin) 
{
  int len=matrix.size();  
  double score=0;
  for(int pos=0, index=begin; pos<len; ++pos,++index)
    score+=matrix[pos]->scoreSingleBase(seq,str,index,seq[index],
					str[index]);

  return score;
}



bool WAM::save(ostream &os)
{
  os.precision(8);
  os << "WAM" << endl;

  os << getSignalType() << endl;
  os << getCutoff() << endl; 
  os << getContextWindowLength() << '\t' 
     << getConsensusOffset() << '\t'
     << getConsensusLength() << endl;  
  os << getStrand() << endl;

  int n=matrix.size(); // n should be equal to ContextWindowLength
  
  for(int i=0 ; i<n ; ++i)
    {
      MarkovChain *chain=matrix[i];
      chain->save(os);
    }
  return true;
}



bool WAM::save(const BOOM::String &filename) 
{
  ofstream os(filename.c_str());
  if(!os.good())
    throw BOOM::String("Error creating file ")+filename+
      "in WAM::save()";
  return save(os);
}



void WAM::load(istream &is)
{
  int contextWindowLength, consensusOffset, consensusLength;
  double cutoff;
  Strand strand;
  SignalType signalType;

  BOOM::String p;
  is >> signalType;
  is >> p; 
  cutoff=p.asDouble();
  is >> contextWindowLength >> consensusOffset >> consensusLength;
  is >> strand;
  
  setSignalType(signalType);
  setStrand(strand);
  setSizes(consensusLength,consensusOffset,contextWindowLength);
  setCutoff(cutoff);

  for(int i=0 ; i< contextWindowLength ; ++i) {
    BOOM::String modelType;
    is >> modelType;
    if(modelType!="MC") 
      throw BOOM::String("Attempt to load an object of type ")+modelType+
	" into a MC";
    matrix.push_back(new MarkovChain(is));
  }
}



SignalSensor *WAM::reverseComplement()
{
  WAM *other=new WAM(getGC(),*this,true);

  int n=matrix.size();
  for(int i=0 ; i<n ; ++i)
    {
      MarkovChain *chain=(MarkovChain*) matrix[i]->reverseComplement();  
      other->matrix.push_front(chain);  // push it in the front
    }

  return other;
}



void WAM::useLogOdds(SignalSensor &nullModel)
{
  WAM &nullWAM=dynamic_cast<WAM&>(nullModel);
  int n=matrix.size();
  for(int i=0 ; i<n ; ++i)
    matrix[i]->useLogOdds(*nullWAM.matrix[i]);
}



void WAM::useLogOdds_anonymous(ContentSensor &nullModel)
{
  int n=matrix.size();
  for(int i=0 ; i<n ; ++i)
    matrix[i]->useLogOdds_anonymous(nullModel);
}
