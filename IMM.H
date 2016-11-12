/****************************************************************
 IMM.H : Interpolated Markov Models
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_IMM_H
#define INCL_IMM_H

#include "ContentSensor.H"
#include "BOOM/Vector.H"
#include "BOOM/StringMap.H"
#include "BOOM/Chi2Table.H"
#include "TrainingSequence.H"

/****************************************************************
 class IMM

 Implements an Interpolated Markov Model as described in: 

    Salzberg SL, Pertea M, Delcher AL, Gardner MJ, Tettelin H: 
    Interpolated Markov models for eukaryotic gene finding.  
    Genomics 59, 1999:24-31.

 There is only one known difference between this implementation and
 that described in the paper: the "minSampleSize" parameter given
 during parameter estimation is used here to enforce a strict lower
 bound on sample size.  Whereas the above paper calls for a rule
 in which probabilities based on sample sizes less than 400 should
 be interpolated with a lower-order estimate, this software does that
 but additionally enforces the rule that any estimate based on a
 sample size less than "minSampleSize" will be completely discarded,
 so that a lower-order estimate must be used. The lower-order estimate
 may then utilize interpolation with a still lower-order model, if
 appropriate.
 ****************************************************************/
class IMM : public ContentSensor
{
  BOOM::Vector<BOOM::StringMap<double>*> *models; // indexed by order
  BOOM::Vector<BOOM::StringMap<int>*> *counts; // indexed by order
  BOOM::Chi2Table chiSquaredTable;
  int N, phase;
  int alphabetSize;
  IMM *revComp;

  double chiTest(BOOM::Vector<int>&,BOOM::Vector<int>&);
  void revCompSeqs(BOOM::Vector<TrainingSequence*>&,
		   BOOM::Vector<TrainingSequence*>&);
  void buildModels(BOOM::Vector<TrainingSequence*> &,int minSampleSize);
  void updateCounts_fw(BOOM::String &,int order,int sequencePhase,
		       int boostCount);
  void updateCounts_rev(BOOM::String &,int order,int seqPhase,
			int boostCount);
  void undefine_fw(BOOM::String &history,BOOM::StringMap<int> &model);
  void undefine_rev(BOOM::String &future,BOOM::StringMap<int> &model);
  void computeProbabilities_fw(int minSampleSize);
  void computeProbabilities_rev(int minSampleSize);
  void load(istream &);
  void interpolate_fw();
  void interpolate_rev();
public:
  IMM(const IMM &);
  IMM(const BOOM::String &filename);
  IMM(istream &,Strand=FORWARD_STRAND);
  IMM(BOOM::Vector<TrainingSequence*> &,int order,int minSampleSize,
      int phase,ContentType,Strand strand=EITHER_STRAND);
  virtual ~IMM();
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
