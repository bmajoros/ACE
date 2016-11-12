/****************************************************************
 SignalSensor.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_SignalSensor_H
#define INCL_SignalSensor_H

#include "BOOM/StringMap.H"
#include "BOOM/String.H"
#include "BOOM/Histogram.H"
#include "SignalType.H"
#include "Signal.H"
#include "BOOM/Strand.H"
class ContentSensor;

class SignalSensor
{
  SignalType signalType; // like ATG, TAG, GT, AG, etc...
  int consensusLength, consensusOffset, contextWindowLength;
  Strand strand;
  double cutoff;
  BOOM::StringMap<char> consensuses;
  GarbageCollector &gc;
  Histogram<double> *probabilityInverter;//logP(seq|sig)->logP(sig|score)
protected:
  SignalSensor(GarbageCollector &);
  SignalSensor(GarbageCollector &,const SignalSensor &,
	       bool reverseComplement=false);
  SignalSensor(const SignalSensor &);
  void setSignalType(SignalType);
  void setSizes(int consensusLength,int consensusOffset,
		int contextWindowLength);
  void setStrand(Strand);
public:
  // Static methods:
  static SignalSensor *load(const BOOM::String &filename,
			    GarbageCollector &);
  static SignalSensor *load(istream &is,GarbageCollector &);

  // Virtual methods:
  virtual SignalPtr detect(const Sequence &,const BOOM::String &,
			 int contextWindowPosition);
  virtual SignalPtr detectWithNoCutoff(const Sequence &,
				       const BOOM::String &,
				       int contextWindowPosition);
  virtual SignalSensor *reverseComplement()=0;
  virtual double getLogP(const Sequence &,const BOOM::String &,int begin)=0;
  virtual bool save(const BOOM::String &filename)=0;
  virtual bool save(ostream &os)=0;
  virtual void addConsensus(const BOOM::String &);
  virtual SignalSensor *clone() const {throw "SignalSensor::clone()";}
  virtual void useLogOdds(SignalSensor &nullModel)=0;
  virtual void useLogOdds_anonymous(ContentSensor &nullModel)=0;

  // Instance methods:
  virtual ~SignalSensor() {}
  GarbageCollector &getGC() {return gc;}
  void useInverseProbabilities(Histogram<double> *);
  void multiplyCutoffBy(double);
  void ignoreCutoff();
  SignalType getSignalType() const;
  virtual bool consensusOccursAt(const BOOM::String &,int index);
  void setCutoff(double);
  double getCutoff() const;
  Strand getStrand() const;
  BOOM::Set<ContentType> &belongsInWhichQueues();
  BOOM::Set<ContentType> &linksBackToWhichQueues();
  SignalPtr getLeftTerminus(double initialTransScore);
  SignalPtr getRightTerminus(int sequenceLength,double finalTransScore);
  int getContextWindowLength() const;
  int getConsensusOffset() const;
  int getConsensusLength() const;
  BOOM::StringMap<char> &getConsensuses() {return consensuses;}
};


#endif
