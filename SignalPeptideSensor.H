/****************************************************************
 SignalPeptideSensor.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalPeptideSensor_H
#define INCL_SignalPeptideSensor_H
#include "BOOM/String.H"
#include "SignalSensor.H"
#include "TrainingSequence.H"
#include "GarbageCollector.H"
#include "CodonTree.H"
#include "WMM.H"

class SignalPeptideSensor : public SignalSensor
{
  WMM *startCodonSensor;
  BOOM::Array1D<CodonTree*> codonModels;
  BOOM::Array1D<int> windowLengths; // only for codonModels

  void revComplementSelf();
  void load(istream &);
  SignalPeptideSensor(GarbageCollector &);
  double scoreCodonModel(CodonTree *,const char *seqStr,int begin,
		       int windowLength);
public:
  SignalPeptideSensor(GarbageCollector &,
		      const SignalPeptideSensor &other,
		      bool reverseComplement=false);
  SignalPeptideSensor(GarbageCollector &,BOOM::String &filename);
  SignalPeptideSensor(GarbageCollector&,istream &);
  SignalPeptideSensor(GarbageCollector &,
		      BOOM::Vector<TrainingSequence*> &,
		      SignalType,
		      int consensusOffset,
		      int consensusLength);
  virtual ~SignalPeptideSensor();
  virtual SignalSensor *reverseComplement();
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  virtual void useLogOdds(SignalSensor &nullModel);
  virtual void useLogOdds_anonymous(ContentSensor &nullModel);
  double getLogP(const Sequence &,const BOOM::String &,int begin);
};


#endif
