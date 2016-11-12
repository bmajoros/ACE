/****************************************************************
 MarkovChain.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_MarkovChain_H
#define INCL_MarkovChain_H

#include "ContentSensor.H"
#include "BOOM/Vector.H"
#include "BOOM/StringMap.H"
#include "TrainingSequence.H"

/*

  NOTE: This Markov chain is not interpolated, but is not fixed-order
        either.  It is a "variable-order" or "back-off" model which
	uses the maximum order possible given the sample size of the
	current history.  This is called a "back off" model in the
	statistical parsing literature.

*/

class MarkovChain : public ContentSensor
{
  BOOM::Vector<BOOM::StringMap<double>*> models; // indexed by order
  int N, phase;
  int alphabetSize;
  MarkovChain *revComp;

  void revCompSeqs(BOOM::Vector<TrainingSequence*>&,
		   BOOM::Vector<TrainingSequence*>&);
  void buildModels(BOOM::Vector<TrainingSequence*> &,int minSampleSize);
  void updateCounts_fw(BOOM::String &,int order,int sequencePhase,
		       int boostCount);
  void updateCounts_rev(BOOM::String &,int order,int seqPhase,
			int boostCount);
  void undefine_fw(BOOM::String &history,BOOM::StringMap<double> &model);
  void undefine_rev(BOOM::String &future,BOOM::StringMap<double> &model);
  void computeProbabilities_fw(int minSampleSize);
  void computeProbabilities_rev(int minSampleSize);
  void load(istream &);
public:
  MarkovChain(const MarkovChain &);
  MarkovChain(const BOOM::String &filename);
  MarkovChain(istream &,Strand=FORWARD_STRAND);
  MarkovChain(BOOM::Vector<TrainingSequence*> &,int order,
	      int minSampleSize,int phase,ContentType,
	      Strand strand=EITHER_STRAND);
  virtual ~MarkovChain();
  virtual double scoreSingleBase(const Sequence &,const BOOM::String &,
				 int index,Symbol,char);
  virtual void scoreSingleBase(const Sequence &,const BOOM::String &,
			       int index,Symbol,char,double &scorePhase0,
			       double &scorePhase1,double &scorePhase2);
  virtual double scoreSubsequence(const Sequence &,const BOOM::String &,
				  int begin,int length,int seqPhase);
  virtual ContentSensor *reverseComplement();
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  virtual void useLogOdds(ContentSensor &nullModel);
  virtual void useLogOdds_anonymous(ContentSensor &nullModel);
  virtual int getOrder() {return N;}
  virtual bool isPhased() {return false;}
  virtual ContentSensor *compile();
};


#endif
