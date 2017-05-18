/****************************************************************
 LogisticSensor.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "LogisticSensor.H"
#include "BOOM/Array1D.H"
#include "BOOM/Exceptions.H"
#include <iostream>
#include <fstream>
#include <math.h>
#include "ContentSensor.H"


const int PRIME_HASH_SIZE=97;
const double log_of_2=log(2.0);
inline float ln(float p) {return log(p)/log_of_2;}



LogisticSensor::LogisticSensor(GarbageCollector &gc,
			       const LogisticSensor &other,bool revComp)
  : matrix(other.matrix), 
    SignalSensor(gc,other,revComp),
    intercept(other.intercept)
{
  // ctor
}



LogisticSensor::LogisticSensor(GarbageCollector &gc)
  : matrix(0,0),
    SignalSensor(gc),
    intercept(0.0)
{
  // protected ctor
}



LogisticSensor::LogisticSensor(GarbageCollector &gc,BOOM::String &filename)
  : matrix(0,0), SignalSensor(gc), intercept(0.0)
{
  // ctor

  ifstream is(filename.c_str());
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="LogisticSensor") 
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      "into a LogisticSensor";
  load(is);
}



LogisticSensor::LogisticSensor(GarbageCollector &gc,istream &is)
  : matrix(0,0), SignalSensor(gc), intercept(0.0)
{
  // ctor

  load(is);
}



LogisticSensor::LogisticSensor(GarbageCollector &gc,
			       BOOM::Vector<TrainingSequence*> &sequences,
			       SignalType signalType,int consensusOffset,
			       int consensusLength)
  : matrix(0,0),
    SignalSensor(gc)
{
  throw "LogisticSensor::LogisticSensor(gc,sequences,signalType,offset,len)";
}



void LogisticSensor::convertToLogs()
{
  throw "LogisticSensor::convertToLogs()";
}



double LogisticSensor::getLogP(const Sequence &seq,const BOOM::String &str,
			       int begin)
{
  int len=matrix.getFirstDim();
  double score=intercept;
  for(int pos=0, index=begin ; pos<len ; ++pos, ++index)
    score+=matrix[pos][seq[index]];
  return score;
}



bool LogisticSensor::save(const BOOM::String &filename) 
{
  ofstream os(filename.c_str());
  if(!os.good())
    throw BOOM::String("Error creating file ")+filename+
      "in LogisticSensor::save()";
  return save(os);
}



bool LogisticSensor::save(ostream &os)
{
  throw "LogisticSensor::save()";
}



void LogisticSensor::load(istream &is)
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
  matrix.setAllTo(0.0);
  String line;
  while(!is.eof()) {
    line.getline(is);
    Vector<String> fields;
    line.getFields(fields);
    if(fields.size()==2) {
      if(fields[0]!="intercept") 
	throw RootException(line+
			" : expecting intercept in LogisticSensor::load()");
      intercept=fields[1].asDouble();
    }
    else if(fields.size()==3) {
      int position=fields[0].asInt();
      int symbol=alphabet.lookup(fields[1][0]);
      matrix[position][symbol]=fields[2].asDouble();
    }
  }  
}



SignalSensor *LogisticSensor::reverseComplement()
{
  LogisticSensor *other=new LogisticSensor(getGC(),*this,true);
  other->revComplementSelf();
  return other;
}



void LogisticSensor::revComplementSelf()
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



void LogisticSensor::swap(char a,char b)
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



void LogisticSensor::useLogOdds(SignalSensor &nullModel)
{
  throw "LogisticSensor::useLogOdds";
}



void LogisticSensor::useLogOdds_anonymous(ContentSensor &nullModel)
{
  throw "ogisticSensor::useLogOdds_anonymous()";
}


float LogisticSensor::divergence(LogisticSensor &other)
{
  throw "LogisticSensor::divergence()";
}
