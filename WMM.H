/****************************************************************
 WMM.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_WMM_H
#define INCL_WMM_H
#include "BOOM/String.H"
#include "BOOM/Array2D.H"
#include "SignalSensor.H"
#include "TrainingSequence.H"
#include "GarbageCollector.H"

class WMM : public SignalSensor
{
protected:
  BOOM::Array2D<float> matrix;

  void load(istream &);
  void revComplementSelf();
  void swap(char,char);
  void convertToLogs();  
  WMM(GarbageCollector &);
public:
  WMM(GarbageCollector &,const WMM &other,bool reverseComplement=false);
  WMM(GarbageCollector &,BOOM::String &filename);
  WMM(GarbageCollector&,istream &);
  
  // This constructor is for training:
  WMM(GarbageCollector &,BOOM::Vector<TrainingSequence*> &,SignalType,
      int consensusOffset,int consensusLength);

  virtual SignalSensor *reverseComplement();
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  virtual void useLogOdds(SignalSensor &nullModel);
  virtual void useLogOdds_anonymous(ContentSensor &nullModel);
  double getLogP(const Sequence &,const BOOM::String &,int begin);
  float divergence(WMM &);
  BOOM::Array2D<float> &getMatrix() {return matrix;}
};


#endif
