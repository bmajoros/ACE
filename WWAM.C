/**************************************************************
 WWAM.C
 mpertea@tigr.org
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include "WWAM.H"
#include <iostream>
#include <fstream>

WWAM::WWAM(GarbageCollector &gc,const WWAM &other,bool revComp)  
  : SignalSensor(gc,other,revComp)
{
  // copy ctor

  if(!revComp) {
    int n=other.matrix.size();  
    for(int i=0;i<n; ++i)
      matrix.push_back(new MarkovChain(*other.matrix[i]));         
  }

  // ### we might want to treat here the revComp==true case
}



WWAM::WWAM(GarbageCollector &gc,BOOM::String &filename)
  : SignalSensor(gc)
{
  // ctor

  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file ")+filename+
		   "in WWAM::WWAM()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="WWAM") 
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      "in into a WWAM";
  load(is);
}



WWAM::WWAM(GarbageCollector &gc,istream &is)
  : SignalSensor(gc)
{
  // ctor

  load(is);
}



WWAM::WWAM(GarbageCollector &gc,BOOM::Vector<TrainingSequence*> &sequences,
	   int order,int minSampleSize, SignalType signalType,
	   int consensusOffset,int consensusLength)
  : SignalSensor(gc)
{
  // ctor

  int samplewindow=2; // ### uses a fixed 5-element averaging window

  BOOM::Vector<TrainingSequence*> examples;

  setStrand(FORWARD_STRAND);
  setSignalType(signalType);
  int windowSize= sequences[0]->getLength() ;
  setSizes(consensusLength,consensusOffset,windowSize); 

  int nosequences=sequences.size();

  for(int i=0; i<windowSize; ++i) {

    int mini=0 > i-samplewindow ? 0 : i-samplewindow;
    int maxi=windowSize-1 < i+samplewindow? windowSize-1 : i+samplewindow;

    for(int j=0;j<nosequences;j++) 
      for(int k=mini;k<=maxi;k++) {
	TrainingSequence *s=new TrainingSequence();
	int min=0 > k-order ? 0 : k-order;
	
	// ### this considers 0-order markov chains for the first base in the window; we might want to try higher orders for this

	sequences[j]->getSubsequence(min,k-min+1,*s);
	examples.push_back(s);
      }
	
    matrix.push_back(new MarkovChain(examples,order,minSampleSize,-1,
				     UNKNOWN_CONTENT_FORWARD));  
    examples.clear();
  }

}



WWAM::~WWAM()
{
  BOOM::Vector<MarkovChain *>::iterator cur=matrix.begin(), end=matrix.end();
  for(; cur!=end ; ++cur)
    delete *cur;
}



double WWAM::getLogP(const Sequence &seq,const BOOM::String &str,int begin) 
{

  int len = matrix.size();  
  double score=0;
  for(int pos=0, index=begin; pos<len; ++pos,++index)
    score+=matrix[pos]->scoreSingleBase(seq,str,index,seq[index],str[index]);

  return score;
}



bool WWAM::save(ostream &os)
{
  os.precision(8);
  os << "WWAM" << endl;

  os << getSignalType() << endl;
  os << getCutoff() << endl; 
  os << getContextWindowLength() << '\t' 
     << getConsensusOffset() << '\t'
     << getConsensusLength() << endl;  
  os << getStrand() << endl;

  int n=matrix.size(); // n should be equal to ContextWindowLength?
  
  for(int i=0 ; i<n ; ++i)
    {
      MarkovChain *chain=matrix[i];
      chain->save(os);
    }
  return true;
}



bool WWAM::save(const BOOM::String &filename) 
{
  ofstream os(filename.c_str());
  if(!os.good())
    throw BOOM::String("Error creating file ")+filename+
      "in WWAM::save()";
  return save(os);
}



void WWAM::load(istream &is)
{
  int contextWindowLength, consensusOffset, consensusLength;
  double cutoff;
  Strand strand;
  SignalType signalType;

  is >> signalType;
  is >> cutoff;
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
	"in into a MC";
    matrix.push_back(new MarkovChain(is));
  }
}



SignalSensor *WWAM::reverseComplement()
{
  WWAM *other=new WWAM(getGC(),*this,true);

  int n=matrix.size();
  for(int i=0 ; i<n ; ++i)
    {
      MarkovChain *chain=(MarkovChain*)matrix[i]->reverseComplement();  
      //other->matrix.push_back(chain);  // push it in the back ->doesn't seem right
      other->matrix.push_front(chain);
    }

  return other;
}



void WWAM::useLogOdds(SignalSensor &nullModel)
{
  WWAM &nullWWAM=dynamic_cast<WWAM&>(nullModel);
  int n=matrix.size();
  for(int i=0 ; i<n ; ++i)
    matrix[i]->useLogOdds(*nullWWAM.matrix[i]);
}



void WWAM::useLogOdds_anonymous(ContentSensor &nullModel)
{
  int n=matrix.size();
  for(int i=0 ; i<n ; ++i)
    matrix[i]->useLogOdds_anonymous(nullModel);
}
