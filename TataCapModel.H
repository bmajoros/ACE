/****************************************************************
 TataCapModel.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TataCapModel_H
#define INCL_TataCapModel_H
#include "SignalSensor.H"
#include "ContentSensor.H"
#include "GarbageCollector.H"
#include "WMM.H"
#include "MarkovChain.H"

class TataCapModel : public SignalSensor
{
  WMM *tataModel, *capModel, *capIntergenicRatioModel;
  MarkovChain *intergenicModel;
  int minSeparation, maxSeparation; // between TATA & CAP context windows
  int tataLength, capLength;

  void revCompSelf();
  TataCapModel();
  void load(istream &);
  void load(const BOOM::String &filename);
public:
  TataCapModel(TataCapModel &other,bool revComp=false);
  TataCapModel(GarbageCollector &gc,const BOOM::String &filename);
  TataCapModel(GarbageCollector &gc,istream &);
  virtual ~TataCapModel();
  virtual bool consensusOccursAt(const BOOM::String &,int index);
  virtual SignalSensor *reverseComplement();
  virtual double getLogP(const Sequence &,const BOOM::String &,int begin);

  // These aren't really needed at the moment:
  virtual bool save(const BOOM::String &) {throw "TataCapModel";}
  virtual bool save(ostream &) {throw "TataCapModel";}
  virtual void addConsensus(const BOOM::String &) {}
  virtual void useLogOdds(SignalSensor &nullModel) 
    {throw "TataCapModel::useLogOdds() not yet implemented";}
  virtual void useLogOdds_anonymous(ContentSensor &nullModel) 
    {throw "TataCapModel::useLogOdds_anonymous() not yet implemented";}
};


#endif
