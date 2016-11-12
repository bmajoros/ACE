/****************************************************************
 WMM.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "WMM.H"
#include "BOOM/Array1D.H"
#include <iostream>
#include <fstream>
#include <math.h>
#include "ContentSensor.H"


const int PRIME_HASH_SIZE=97;
const double log_of_2=log(2.0);
inline float ln(float p) {return log(p)/log_of_2;}



WMM::WMM(GarbageCollector &gc,const WMM &other,bool revComp)
  : matrix(other.matrix), 
    SignalSensor(gc,other,revComp)
{
  // ctor
}



WMM::WMM(GarbageCollector &gc)
  : matrix(0,0),
    SignalSensor(gc)
{
  // protected ctor
}



WMM::WMM(GarbageCollector &gc,BOOM::String &filename)
  : matrix(0,0),
    SignalSensor(gc)
{
  // ctor

  ifstream is(filename.c_str());
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="WMM") 
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      "into a WMM";
  load(is);
}



WMM::WMM(GarbageCollector &gc,istream &is)
  : matrix(0,0),
    SignalSensor(gc)
{
  // ctor

  load(is);
}



WMM::WMM(GarbageCollector &gc,BOOM::Vector<TrainingSequence*> &sequences,
	 SignalType signalType,int consensusOffset,int consensusLength)
  : matrix(0,0),
    SignalSensor(gc)
{
  /*
    This constructor performs training of the WMM
   */

  setStrand(FORWARD_STRAND);
  setSignalType(signalType);

  float pseudocount=0;
  int n=sequences.size();
  int len=sequences[0]->getLength();
  setSizes(consensusLength,consensusOffset,len);
  int nAlpha=alphabet.getNumElements();
  matrix.resize(len,nAlpha);

  matrix.setAllTo(pseudocount);

  BOOM::FloatArray1D effectiveSize(len);
  effectiveSize.setAllTo(0);
  for(int i=0 ; i<n ; ++i)
    {
      TrainingSequence &seq=*sequences[i];
      int l=seq.getLength();
      if(l!=len) throw BOOM::String("length mismatch in WMM: ")+len+" vs "+l;
      for(int pos=0 ; pos<len ; ++pos)
	{
	  Symbol s=seq[pos];
	  int count=seq.getBoostCount();
	  matrix[pos][s]+=count;
	  effectiveSize[pos]+=count;
	}
    }
  for(int pos=0 ; pos<len ; ++pos)
    for(Symbol s=0 ; s<nAlpha ; ++s)
      matrix[pos][s]/=effectiveSize[pos];

  convertToLogs();  
}



void WMM::convertToLogs()
{
  int maxX=matrix.getFirstDim(), maxY=matrix.getSecondDim();
  for(int x=0 ; x<maxX ; ++x)
    for(int y=0 ; y<maxY ; ++y)
      matrix[x][y]=log(matrix[x][y]);
}



double WMM::getLogP(const Sequence &seq,const BOOM::String &str,int begin)
{
  int len=matrix.getFirstDim();
  double score=0;
  for(int pos=0, index=begin ; pos<len ; ++pos, ++index)
    score+=matrix[pos][seq[index]];
  return score;
}



bool WMM::save(const BOOM::String &filename) 
{
  ofstream os(filename.c_str());
  if(!os.good())
    throw BOOM::String("Error creating file ")+filename+
      "in WMM::save()";
  return save(os);
}



bool WMM::save(ostream &os)
{
  os.precision(8);
  int maxX=matrix.getFirstDim(), maxY=matrix.getSecondDim();

  os << "WMM" << endl;
  os << getSignalType() << endl;
  os << getCutoff() << "\t" << maxX << "\t" << maxY << endl;
  os << getContextWindowLength() << '\t' 
     << getConsensusOffset() << '\t'
     << getConsensusLength() << endl;
  os << getStrand() << endl;

  for(int x=0 ; x<maxX ; ++x)
    {
      for(int y=0 ; y<maxY ; ++y)
	os << matrix[x][y] << "\t";
      os << endl;
    }
  return true;
}



void WMM::load(istream &is)
{
  int maxX, maxY, contextWindowLength, consensusOffset, consensusLength;
  double cutoff;
  Strand strand;
  SignalType signalType;

  BOOM::String p;
  is >> signalType;
  is >> p >> maxX >> maxY;
  cutoff=p.asDouble();
  is >> contextWindowLength >> consensusOffset >> consensusLength;
  is >> strand;
  
  setSignalType(signalType);
  setStrand(strand);
  setSizes(consensusLength,consensusOffset,contextWindowLength);
  setCutoff(cutoff);

  matrix.resize(maxX,maxY);
  for(int x=0 ; x<maxX ; ++x)
    for(int y=0 ; y<maxY ; ++y)
      {
	is >> p;
	matrix[x][y]=p.asDouble();
      }  
}



SignalSensor *WMM::reverseComplement()
{
  WMM *other=new WMM(getGC(),*this,true);
  other->revComplementSelf();
  return other;
}



void WMM::revComplementSelf()
{
  // Complement:
  swap('A','T');
  swap('C','G');

  // Reverse:
  int seqLen=matrix.getFirstDim();
  int nAlpha=alphabet.getNumElements();
  int halfLen=seqLen/2;
  for(Symbol s=0 ; s<nAlpha ; ++s)
    for(int i=0 ; i<halfLen ; ++i)
      {
	float temp=matrix[i][s];
	matrix[i][s]=matrix[seqLen-i-1][s];
	matrix[seqLen-i-1][s]=temp;
      }
}



void WMM::swap(char a,char b)
{
  int ia=alphabet.lookup(a), ib=alphabet.lookup(b);
  int seqLen=matrix.getFirstDim();
  for(int i=0 ; i<seqLen ; ++i)
    {
      float temp=matrix[i][ia];
      matrix[i][ia]=matrix[i][ib];
      matrix[i][ib]=temp;
    }
}



void WMM::useLogOdds(SignalSensor &nullModel)
{
  WMM &nullWMM=dynamic_cast<WMM&>(nullModel);
  int seqLen=matrix.getFirstDim();
  int nAlpha=matrix.getSecondDim();
  for(int i=0 ; i<seqLen ; ++i)
    for(Symbol s=0 ; s<nAlpha ; ++s)
      matrix[i][s]-=nullWMM.matrix[i][s];
}



void WMM::useLogOdds_anonymous(ContentSensor &nullModel)
{
  int seqLen=matrix.getFirstDim();
  int nAlpha=matrix.getSecondDim();
  for(int i=0 ; i<seqLen ; ++i)
    for(Symbol s=0 ; s<nAlpha ; ++s)
      {
	char str[2];
	str[0]=alphabet.lookup(s);
	str[1]='\0';
	Sequence seq(str,alphabet);
	matrix[i][s]-=nullModel.scoreSingleBase(seq,str,0,s,str[0]);
      }
}


float WMM::divergence(WMM &other)
{
  float div=0;
  int maxX=matrix.getFirstDim(), maxY=matrix.getSecondDim();
  for(int x=0 ; x<maxX ; ++x)
    for(int y=0 ; y<maxY ; ++y)
      {
	float p=matrix[x][y], q=matrix[x][y];
	/*if(q)*/ div+=p*fabs(ln(p/q));
      }
  return div;
}
