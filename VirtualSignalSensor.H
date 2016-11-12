/****************************************************************
 VirtualSignalSensor.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VirtualSignalSensor_H
#define INCL_VirtualSignalSensor_H
#include <iostream>
#include "SignalSensor.H"
using namespace std;
using namespace BOOM;


/*
  Allows creation of virtual signals such as LEFT_TERMINUS and
  RIGHT_TERMINUS, as well as signals created from RNA-seq or
  annotations, such as TSS (transcription start site) or TES
  (transcript end site).
 */
class VirtualSignalSensor : public SignalSensor {
public:
  VirtualSignalSensor(SignalType,GarbageCollector &);
  virtual bool consensusOccursAt(const BOOM::String &,int index);
  virtual SignalPtr detect(const Sequence &,const BOOM::String &,
			 int contextWindowPosition);
  virtual SignalPtr detectWithNoCutoff(const Sequence &,
				       const BOOM::String &,
				       int contextWindowPosition);
  virtual SignalSensor *reverseComplement();
  virtual double getLogP(const Sequence &,const BOOM::String &,int begin);
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  virtual void useLogOdds(SignalSensor &nullModel);
  virtual void useLogOdds_anonymous(ContentSensor &nullModel);
};

#endif

