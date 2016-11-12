/****************************************************************
 Fast3PMC.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_Fast3PMC_H
#define INCL_Fast3PMC_H

#include "BOOM/Alphabet.H"
#include "BOOM/File.H"
#include "ContentSensor.H"
#include "FastMarkovChain.H"
#include "ThreePeriodicMarkovChain.H"
#include "ThreePeriodicIMM.H"

class Fast3PMC : public ContentSensor
{
  FastMarkovChain *chains[3];
  Fast3PMC *revComp;

  void load(const BOOM::String &filename);
  void load(BOOM::File &);
  void compileFrom(ThreePeriodicMarkovChain &);
  void compileFrom(ThreePeriodicIMM &);
  Fast3PMC(ContentType);
public:
  Fast3PMC(ThreePeriodicMarkovChain &);
  Fast3PMC(ThreePeriodicIMM &);
  Fast3PMC(const BOOM::String &filename);
  Fast3PMC(BOOM::File &);
  virtual int getOrder() {return chains[0]->getOrder();}
  virtual bool isPhased() {return true;}
  virtual double scoreSingleBase(const Sequence &,const BOOM::String &,
				 int index,Symbol,char);
  virtual void scoreSingleBase(const Sequence &,const BOOM::String &,
			       int index,Symbol,char,double &scorePhase0,
			       double &scorePhase1,double &scorePhase2);
  // NOTE: before calling scoreSubsequence() you should call reset() to
  // place the model into its initial state...unless you really
  // know what you're doing...
  virtual double scoreSubsequence(const Sequence &,const BOOM::String &,
				  int begin,int length,int phase);
  virtual ContentSensor *reverseComplement();
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  bool save(BOOM::File &);
  virtual void useLogOdds(ContentSensor &nullModel) 
    {throw "Fast3PMC::useLogOdds() not yet implemented";}
  virtual void useLogOdds_anonymous(ContentSensor &nullModel) 
    {throw "Fast3PMC::useLogOdds_anonymous() not yet implemented";}
  virtual void reset(const Sequence &,const BOOM::String &,int pos);
};


#endif
