/****************************************************************
 LogisticSensor.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_LogisticSensor_H
#define INCL_LogisticSensor_H
#include "BOOM/String.H"
#include "BOOM/Array2D.H"
#include "SignalSensor.H"
#include "TrainingSequence.H"
#include "GarbageCollector.H"

class LogisticSensor : public SignalSensor
{
protected:
  BOOM::Array2D<float> matrix;
  BOOM::Array2D<float> Pmatrix;
  float intercept;
  void load(istream &);
  void revComplementSelf();
  void swap(char,char);
  void convertToLogs();  
  LogisticSensor(GarbageCollector &);
public:
  LogisticSensor(GarbageCollector &,const LogisticSensor &other,
		 bool reverseComplement=false);
  LogisticSensor(GarbageCollector &,BOOM::String &filename);
  LogisticSensor(GarbageCollector&,istream &);
  
  // This constructor is for training:
  LogisticSensor(GarbageCollector &,BOOM::Vector<TrainingSequence*> &,
		 SignalType,int consensusOffset,int consensusLength);

  virtual SignalSensor *reverseComplement();
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  virtual void useLogOdds(SignalSensor &nullModel);
  virtual void useLogOdds_anonymous(ContentSensor &nullModel);
  double getLogP(const Sequence &,const BOOM::String &,int begin);
  double getRawScore(const Sequence &,const BOOM::String &,int begin);
  float divergence(LogisticSensor &);
  BOOM::Array2D<float> &getMatrix() {return matrix;}

  // Specific to this class:
  double applyLogistic(double BX); // P(splice site | sequence)
};


#endif
