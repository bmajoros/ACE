/****************************************************************
 ThreePeriodicMarkovChain.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_ThreePeriodicMarkovChain_H
#define INCL_ThreePeriodicMarkovChain_H

#include "ContentSensor.H"
#include "MarkovChain.H"
class Fast3PMC;

class ThreePeriodicMarkovChain : public ContentSensor
{
  friend class Fast3PMC;
  MarkovChain *chains[3];
  
  void buildModels(BOOM::Vector<TrainingSequence*> &,int minSampleSize,
		   int order);
  void load(istream &);
public:
  ThreePeriodicMarkovChain(const BOOM::String &filename);
  ThreePeriodicMarkovChain(istream &);
  ThreePeriodicMarkovChain(Strand,ContentType);
  ThreePeriodicMarkovChain(BOOM::Vector<TrainingSequence*> &,int order,
			   int minSampleSize,ContentType);
  virtual ~ThreePeriodicMarkovChain();
  virtual int getOrder() {return chains[0]->getOrder();}
  virtual bool isPhased() {return true;}
  virtual double scoreSingleBase(const Sequence &,const BOOM::String &,
				 int index,Symbol,char);
  virtual void scoreSingleBase(const Sequence &,const BOOM::String &,
			       int index,Symbol,char,double &scorePhase0,
			       double &scorePhase1,double &scorePhase2);
  virtual double scoreSubsequence(const Sequence &,const BOOM::String &,
				  int begin,int length,int seqPhase);
  virtual ContentSensor *reverseComplement();
  virtual bool save(const BOOM::String &filename);
  virtual void useLogOdds(ContentSensor &nullModel);
  virtual void useLogOdds_anonymous(ContentSensor &nullModel);
  virtual bool save(ostream &os);
  virtual ContentSensor *compile();
};


#endif
