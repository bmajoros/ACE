/****************************************************************
 ThreePeriodicIMM.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_ThreePeriodicIMM_H
#define INCL_ThreePeriodicIMM_H

#include "ContentSensor.H"
#include "IMM.H"
class Fast3PMC;

class ThreePeriodicIMM : public ContentSensor
{
  friend class Fast3PMC;
  IMM *chains[3];
  ThreePeriodicIMM *revComp;
  
  void buildModels(BOOM::Vector<TrainingSequence*> &,int minSampleSize,
		   int order);
  void load(istream &);
public:
  ThreePeriodicIMM(const BOOM::String &filename);
  ThreePeriodicIMM(istream &);
  ThreePeriodicIMM(Strand,ContentType);
  ThreePeriodicIMM(BOOM::Vector<TrainingSequence*> &,int order,
			   int minSampleSize,ContentType);
  virtual ~ThreePeriodicIMM();
  IMM *&getChain(int i) { return chains[i]; }
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
  virtual int getOrder() {return chains[0]->getOrder();}
  virtual bool isPhased() {return true;}
  virtual ContentSensor *compile();
};


#endif
