/****************************************************************
 BranchAcceptor.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_BranchAcceptor_H
#define INCL_BranchAcceptor_H
#include "SignalSensor.H"
#include "TrainingSequence.H"
#include "WWAM.H"
#include "WAM.H"

class BranchAcceptor : public SignalSensor
{
                     // GENSCAN uses:
  WAM *acceptor;     //   15 bp around acceptor consensus (AG)
  WWAM *branchPoint; //   19 bp immediately upstream from acceptor WAM

  void load(istream &);
public:
  BranchAcceptor(GarbageCollector &,const BranchAcceptor &, 
		 bool reverseComplement=false);
  BranchAcceptor(GarbageCollector &,BOOM::String &filename);
  BranchAcceptor(GarbageCollector &,istream &);
  BranchAcceptor(GarbageCollector &,
		 BOOM::Vector<TrainingSequence*> &,
		 int branchPointOrder,
		 int acceptorOrder,
		 int branchContextLength, // length of WWAM
		 int minSampleSize,
		 int consensusOffset,     // from beginning of whole sensor
		 SignalType);             // AG or NEG_AG
  virtual ~BranchAcceptor();
  
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  double getLogP(const Sequence &,const BOOM::String &,int begin);
  virtual SignalSensor *reverseComplement();
  virtual void useLogOdds(SignalSensor &nullModel);
  virtual void useLogOdds_anonymous(ContentSensor &nullModel);
};


#endif
