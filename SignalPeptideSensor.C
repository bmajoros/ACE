/****************************************************************
 SignalPeptideSensor.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "SignalPeptideSensor.H"
#include <iostream>


SignalPeptideSensor::SignalPeptideSensor(GarbageCollector &gc)
  : SignalSensor(gc),
    startCodonSensor(NULL)
{
  // ctor
}



SignalPeptideSensor::SignalPeptideSensor(GarbageCollector &gc,
					 const SignalPeptideSensor &other,
					 bool revComp)
  : SignalSensor(gc,other,revComp),
    startCodonSensor(new WMM(gc,*other.startCodonSensor,revComp)),
    windowLengths(other.windowLengths.size()),
    codonModels(other.codonModels.size())
{
  // copy ctor

  int n=other.windowLengths.size();
  for(int i=0 ; i<n ; ++i)
    {
      windowLengths[i]=other.windowLengths[i];
      codonModels[i]=other.codonModels[i]->clone();
    }
}



SignalPeptideSensor::SignalPeptideSensor(GarbageCollector &gc,
					 BOOM::String &filename)
  : SignalSensor(gc)
{
  // ctor

  ifstream is(filename.c_str());
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="SignalPeptide") 
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      "into a SignalPeptide sensor";
  load(is);
}



SignalPeptideSensor::SignalPeptideSensor(GarbageCollector &gc,
					 istream &is)
  : SignalSensor(gc)
{
  // ctor

  load(is);
}



SignalPeptideSensor::SignalPeptideSensor(GarbageCollector &gc,
					 BOOM::Vector<TrainingSequence*> &T,
					 SignalType signalType,
					 int consensusOffset,
					 int consensusLength)
  : SignalSensor(gc)
{
  // training ctor

  throw "Not implemented here -- use train-signal-peptide-model.C";
}



SignalPeptideSensor::~SignalPeptideSensor()
{
  // dtor

  delete startCodonSensor;
  int n=codonModels.size();
  for(int i=0 ; i<n ; ++i)
    delete codonModels[i];
}



SignalSensor *SignalPeptideSensor::reverseComplement()
{
  SignalPeptideSensor *other=new SignalPeptideSensor(getGC(),*this,true);
  other->revComplementSelf();
  return other;
}



bool SignalPeptideSensor::save(const BOOM::String &filename)
{
  throw "Not implemented here -- use train-signal-peptide-model.C";
}



bool SignalPeptideSensor::save(ostream &os)
{
  throw "Not implemented here -- use train-signal-peptide-model.C";
}



void SignalPeptideSensor::useLogOdds(SignalSensor &nullModel)
{
  throw "Not implemented";
}



void SignalPeptideSensor::useLogOdds_anonymous(ContentSensor &nullModel)
{
  throw "Not implemented";
}



double SignalPeptideSensor::getLogP(const Sequence &seq,
				    const BOOM::String &str,
				    int begin)
{
  switch(getStrand())
    {
    case FORWARD_STRAND:
      {
	double logP=startCodonSensor->getLogP(seq,str,begin);
	int pos=begin+startCodonSensor->getContextWindowLength();
	int n=codonModels.size();
	const char *charPtr=str.c_str();
	for(int i=0 ; i<n ; ++i)
	  {
	    CodonTree *codonModel=codonModels[i];
	    int windowLength=windowLengths[i];
	    logP+=scoreCodonModel(codonModel,charPtr,pos,windowLength);
	    pos+=windowLength;
	  }
	return logP;
      }
      break;
    case REVERSE_STRAND:
      {
	int pos=begin+getContextWindowLength()-
	  startCodonSensor->getContextWindowLength();
	double logP=startCodonSensor->getLogP(seq,str,pos);
	int n=codonModels.size();
	const char *charPtr=str.c_str();
	for(int i=0 ; i<n ; ++i)
	  {
	    CodonTree *codonModel=codonModels[i];
	    int windowLength=windowLengths[i];
	    logP+=scoreCodonModel(codonModel,charPtr,pos,windowLength);
	    pos-=windowLength;
	  }
	return logP;
      }
      break;
    }
}



double SignalPeptideSensor::scoreCodonModel(CodonTree *codonModel,
					    const char *seqStr,
					    int begin,
					    int windowLength)
{
  double logP=0;
  const char *p=seqStr+begin;
  for(int i=0 ; i<windowLength ; ++i, p+=3)
    logP+=codonModel->scoreCodon(p);
  return logP;
}



void SignalPeptideSensor::revComplementSelf()
{
  startCodonSensor=(WMM*) startCodonSensor->reverseComplement();
  int n=codonModels.size();
  for(int i=0 ; i<n ; ++i)
    codonModels[i]->revComplementSelf();
}



void SignalPeptideSensor::load(istream &is)
{
  BOOM::String wmm;
  is>>wmm;
  if(wmm!="WMM") throw "start codon model must be WMM when used within a SignalPeptideSensor";
  startCodonSensor=new WMM(getGC(),is);

  setSignalType(ATG);
  setStrand(FORWARD_STRAND);
  ignoreCutoff(); // ### 3/17/05 WHM

  int contextWindowLength=startCodonSensor->getContextWindowLength();
  int numFields;
  is>>numFields;
  codonModels.resize(numFields);
  windowLengths.resize(numFields);
  for(int i=0 ; i<numFields ; ++i)
    {
      is>>windowLengths[i];
      codonModels[i]=new CodonTree(is);
      contextWindowLength+=3*windowLengths[i];
    }

  setSizes(3,startCodonSensor->getConsensusOffset(),
	   contextWindowLength);
}


