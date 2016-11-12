/*******************************************************************
 WAM.H : Weight Array Matrix
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ******************************************************************/
#ifndef INCL_WAM_H
#define INCL_WAM_H
#include "BOOM/String.H"
#include "BOOM/Vector.H"
#include "MarkovChain.H"
#include "SignalSensor.H"
#include "TrainingSequence.H"
#include "GarbageCollector.H"

class WAM : public SignalSensor
{
  BOOM::Vector<MarkovChain *> matrix;

  void load(istream &);
public:
  WAM(GarbageCollector &,const WAM &, bool reverseComplement=false);
  WAM(GarbageCollector &,BOOM::String &filename);
  WAM(GarbageCollector &,istream &);
  WAM(GarbageCollector &,BOOM::Vector<TrainingSequence*> &, int order, 
      int minSampleSize,SignalType,int consensusOffset,int consensusLength); 
  virtual ~WAM();
  
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  double getLogP(const Sequence &,const BOOM::String &,int begin);
  virtual void useLogOdds(SignalSensor &nullModel);
  virtual void useLogOdds_anonymous(ContentSensor &nullModel);
  virtual SignalSensor *reverseComplement();
};

#endif
