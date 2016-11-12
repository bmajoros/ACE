/****************************************************************
 ContentSensor.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "ContentSensor.H"
#include <iostream>
#include <fstream>
#include "MarkovChain.H"
#include "ThreePeriodicMarkovChain.H"
#include "FastMarkovChain.H"
#include "Fast3PMC.H"
#include "IMM.H"
#include "ThreePeriodicIMM.H"
#include "SignalQueue.H"
#include "BOOM/PureDnaAlphabet.H"


BOOM::Regex ContentSensor::binmodRegex("binmod$");
const int ContentSensor::NO_PHASE=-1;


ContentSensor::ContentSensor()
{
}



void ContentSensor::setContentType(ContentType t)
{
  contentType=t;
}



void ContentSensor::setStrand(Strand s)
{
  strand=s;
}



ContentSensor *ContentSensor::loadBinary(const BOOM::String &filename)
{
  BOOM::File file(filename);
  BOOM::String modelType;
  file >> modelType;
  file.rewind();

  if(modelType=="FastMC") return new FastMarkovChain(file);
  if(modelType=="Fast3PMC") return new Fast3PMC(file);

  throw BOOM::String("Unrecognized model type in ContentSensor::loadBinary(): ")+
    modelType;
}



ContentSensor *ContentSensor::load(const BOOM::String &filename)
{
  if(binmodRegex.search(filename))
    return loadBinary(filename);

  ifstream is(filename.c_str());
  if(!is.good())
    throw BOOM::String("Error opening file ")+filename+
      " in ContentSensor::load()";
  return load(is);    
}



ContentSensor *ContentSensor::load(istream &is)
{
  BOOM::String modelType;
  is >> modelType;

  if(modelType=="MC")
    return new MarkovChain(is);
  if(modelType=="3P")
    return new ThreePeriodicMarkovChain(is);
  if(modelType=="IMM")
    return new IMM(is);
  if(modelType=="3PIMM")
    return new ThreePeriodicIMM(is);

  /*
  if(modelType=="HMM")
    return new HiddenMarkovModel(is,alphabet);
  if(modelType=="NonstationaryMarkovChain")
    return new NonstationaryMarkovChain(is,alphabet);
  if(modelType=="CodingPotential")
    return new CodingPotential(is,alphabet);
  */

  throw BOOM::String("Unrecognized model type in ContentSensor::load(): ")+
    modelType;
}



ContentType ContentSensor::getContentType() const
{
  return contentType;
}



Strand ContentSensor::getStrand() const
{
  return strand;
}



bool ContentSensor::isCoding() const
{
  return ::isCoding(contentType);
}



bool ContentSensor::isIntron() const
{
  return ::isIntron(contentType);
}



bool ContentSensor::isIntergenic() const
{
  return ::isIntergenic(contentType);
}



bool ContentSensor::isUTR() const
{
  return ::isUTR(contentType);
}


void ContentSensor::addQueue(SignalQueue *q)
{
  signalQueues.insert(q);
}



BOOM::Set<SignalQueue*> &ContentSensor::getSignalQueues()
{
  return signalQueues;
}



void ContentSensor::computeScores(const Sequence &seq,Array1D<double> &scores)
{
  String str=seq(PureDnaAlphabet::global());
  int L=str.length();
  scores.resize(L);
  for(int i=0 ; i<L ; ++i)
    scores[i]=scoreSingleBase(seq,str,i,seq[i],str[i]);
}



void ContentSensor::computeScores(const Sequence &seq,Array1D<double> &phase0,
				  Array1D<double> &phase1,
				  Array1D<double> &phase2)
{
  String str=seq(PureDnaAlphabet::global());
  int L=str.length();
  phase0.resize(L); phase1.resize(L); phase2.resize(L);
  for(int i=0 ; i<L ; ++i)
    scoreSingleBase(seq,str,i,seq[i],str[i],phase0[i],phase1[i],phase2[i]);
}



