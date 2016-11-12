/****************************************************************
 FastMarkovChain.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_FastMarkovChain_H
#define INCL_FastMarkovChain_H

#include "BOOM/Array2D.H"
#include "BOOM/Array1D.H"
#include "BOOM/Alphabet.H"
#include "BOOM/File.H"
#include "ContentSensor.H"
class MarkovChainCompiler;

class FastMarkovChain : public ContentSensor
{
  friend class MarkovChainCompiler;
  BOOM::Array2D<int> transitionMatrix;
  BOOM::Array1D<double> probabilityVector;
  BOOM::Array1D<int> statesPerOrder;
  int numStates;
  int order;
  int alphabetSize;
  int state;
  int seqLen;
  FastMarkovChain *revComp;

  void computeNumberOfStates();
  void load(const BOOM::String &filename);
  void load(BOOM::File &);
  void loadNonRecurs(BOOM::File &);
  int stringToState(const BOOM::String &s);
  BOOM::String stateToString(int state);
  bool saveNonRecurs(BOOM::File &);
public:
  FastMarkovChain(int order,ContentType);
  FastMarkovChain(const BOOM::String &filename);
  FastMarkovChain(BOOM::File &);
  virtual int getOrder() {return order;}
  virtual bool isPhased() {return false;}
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os); // don't use this
  virtual bool save(BOOM::File &);
  virtual double scoreSingleBase(const Sequence &,const BOOM::String &,
				 int index,Symbol,char);
  virtual void scoreSingleBase(const Sequence &,const BOOM::String &,
			       int index,Symbol,char,double &scorePhase0,
			       double &scorePhase1,
			       double &scorePhase2);

  // NOTE: before calling scoreSubsequence() you should call reset() to
  // place the model into its initial state...unless you really
  // know what you're doing...
  virtual double scoreSubsequence(const Sequence &,const BOOM::String &,
				  int begin,int length,int phase);

  virtual ContentSensor *reverseComplement();
  void setTrans(int state1,Symbol s,int state2);
  void setProb(int state,double prob);
  int statesOfOrder(int order);
  virtual void reset(const Sequence &,const BOOM::String &,int pos);
  BOOM::String getState(); // useful for debugging
  void gotoState(const BOOM::String &);
  virtual void useLogOdds(ContentSensor &nullModel) 
    {throw "FastMarkovChain::useLogOdds() not implemented yet";}
  virtual void useLogOdds_anonymous(ContentSensor &nullModel) 
    {throw "FastMarkovChain::useLogOdds_anonymous() not implemented yet";}
};


#endif
