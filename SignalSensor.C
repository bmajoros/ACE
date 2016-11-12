/****************************************************************
 SignalSensor.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "SignalSensor.H"
#include <iostream>
#include <fstream>
#include "BOOM/ProteinTrans.H"
#include "BOOM/Constants.H"
#include "WMM.H"
#include "SignalTypeProperties.H"
#include "WAM.H"
#include "WWAM.H"
#include "MddTree.H"
#include "BranchAcceptor.H"
#include "SignalPeptideSensor.H"

const int PRIME_HASH_SIZE=97;


SignalSensor::SignalSensor(const SignalSensor &other)
  : signalType(other.signalType),
    consensusLength(other.consensusLength),
    consensusOffset(other.consensusOffset),
    contextWindowLength(other.contextWindowLength),
    strand(other.strand),
    cutoff(other.cutoff),
    gc(other.gc),
    consensuses(other.consensuses),
    probabilityInverter(NULL)
{
}



SignalSensor::SignalSensor(GarbageCollector &gc)
  : consensusLength(0),
    consensuses(PRIME_HASH_SIZE),
    cutoff(0),
    gc(gc),
    probabilityInverter(NULL)
{
}



SignalSensor::SignalSensor(GarbageCollector &gc,const SignalSensor &other,
			   bool revComp)
  : consensuses(PRIME_HASH_SIZE),
    gc(gc),
    probabilityInverter(NULL)
{
  contextWindowLength=other.contextWindowLength;
  consensusLength=other.consensusLength;
  cutoff=other.cutoff;
  
  if(revComp)
    {
      consensusOffset=
	contextWindowLength-other.consensusOffset-consensusLength;
      strand=complement(other.strand);
      signalType=::reverseComplement(other.signalType);
      BOOM::StringMap<char>::const_iterator cur=other.consensuses.begin(), 
	end=other.consensuses.end();
      for(; cur!=end ; ++cur)
	{
	  const StringMapElem<char> &elem=*cur;
	  BOOM::String reverseConsensus=
	    BOOM::ProteinTrans::reverseComplement(elem.first);
	  consensuses.lookup(reverseConsensus.c_str(),elem.len)=char(1);
	}
    }
  else
    {
      consensusOffset=other.consensusOffset;
      strand=other.strand;
      signalType=other.signalType;
      consensuses=other.consensuses;
    }
}



void SignalSensor::multiplyCutoffBy(double m)
{
  cutoff*=m;
}



SignalSensor *SignalSensor::load(const BOOM::String &filename,
				 GarbageCollector &gc)
{
  ifstream is(filename.c_str());
  if(!is.good())
    throw BOOM::String("Error opening file ")+filename+
      " in SignalSensor::load()";
  return load(is,gc);
}



SignalSensor *SignalSensor::load(istream &is,GarbageCollector &gc)
{
  BOOM::String modelType;
  is >> modelType;

  if(modelType=="WMM")
    return new WMM(gc,is);
  if(modelType=="WAM")
    return new WAM(gc,is);
  if(modelType=="WWAM")
    return new WWAM(gc,is);
  if(modelType=="MDD")
    return new MddTree(is,gc);
  if(modelType=="BranchAcceptor")
    return new BranchAcceptor(gc,is);
  if(modelType=="SignalPeptide")
    return new SignalPeptideSensor(gc,is);

  throw BOOM::String("Unrecognized model type in SignalSensor::load(): ")+
    modelType;
}



SignalType SignalSensor::getSignalType() const
{
  return signalType;
}



void SignalSensor::addConsensus(const BOOM::String &s)
{
  int len=s.length();
  if(consensusLength>0 && consensusLength!=len)
    throw BOOM::String(
       "Consensus lengths differ in SignalSensor::addConsensus");
  
  consensusLength=len;
  consensuses.lookup(s.c_str(),len)=char(1);
}



bool SignalSensor::consensusOccursAt(const BOOM::String &str,int index)
{
  return consensuses.isDefined(str.c_str(),index,consensusLength);
}



void SignalSensor::setSizes(int consensusLength,int consensusOffset,
			    int contextWindowLength)
{
  this->consensusLength=consensusLength;
  this->consensusOffset=consensusOffset;
  this->contextWindowLength=contextWindowLength;
}



double SignalSensor::getCutoff() const
{
  return cutoff;
}



void SignalSensor::setStrand(Strand s)
{
  strand=s;
}



Strand SignalSensor::getStrand() const
{
  return strand;
}



void SignalSensor::setSignalType(SignalType t)
{
  signalType=t;
}



int SignalSensor::getContextWindowLength() const
{
  return contextWindowLength;
}



int SignalSensor::getConsensusOffset() const
{
  return consensusOffset;
}



int SignalSensor::getConsensusLength() const
{
  return consensusLength;
}



void SignalSensor::setCutoff(double c)
{
  cutoff=c;
}



BOOM::Set<ContentType> &SignalSensor::belongsInWhichQueues()
{
  return SignalTypeProperties::global.belongsInWhichQueues(signalType);
}



BOOM::Set<ContentType> &SignalSensor::linksBackToWhichQueues()
{
  return SignalTypeProperties::global.linksBackToWhichQueues(signalType);
}



SignalPtr SignalSensor::getLeftTerminus(double initialTransScore)
{
  SignalPtr dummy=
    new Signal(-getContextWindowLength(),initialTransScore,*this,gc,
	       signalType);
  return dummy;
}



SignalPtr SignalSensor::getRightTerminus(int sequenceLength,
					 double finalTransScore)
{
  SignalPtr dummy=
    new Signal(sequenceLength,finalTransScore,*this,gc,signalType);
  return dummy;
}



SignalPtr SignalSensor::detect(const Sequence &seq,const BOOM::String &str,
			     int contextWindowPosition)
{
  // First, check whether a legal consensus sequence occurs at the
  // appropriate offset
  int consensusPosition=contextWindowPosition+consensusOffset;
  if(!consensusOccursAt(str,consensusPosition)) return NULL;

  double score=getLogP(seq,str,contextWindowPosition);
  if(score>=cutoff)
    {
      //if(probabilityInverter) score=probabilityInverter->getP(score);
      return new Signal(contextWindowPosition,score,*this,gc,
			signalType);
    }
  return NULL;
}



SignalPtr SignalSensor::detectWithNoCutoff(const Sequence &seq,
					   const BOOM::String &str,
					   int contextWindowPosition)
{
  int consensusPosition=contextWindowPosition+consensusOffset;
  if(!consensusOccursAt(str,consensusPosition)) return NULL;
  double score=getLogP(seq,str,contextWindowPosition);
  //if(probabilityInverter) score=probabilityInverter->getP(score);
  return new Signal(contextWindowPosition,score,*this,gc,signalType);
}



void SignalSensor::ignoreCutoff()
{
  cutoff=NEGATIVE_INFINITY;
}



void SignalSensor::useInverseProbabilities(Histogram<double> *h)
{
  probabilityInverter=h;
}


