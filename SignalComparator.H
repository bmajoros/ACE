/****************************************************************
 SignalComparator.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalComparator_H
#define INCL_SignalComparator_H
#include "BOOM/Comparator.H"
#include "Signal.H"
#include "DiscreteDistribution.H"
#include "ContentType.H"


class IgnorantComparator : public BOOM::Comparator<SignalPtr>
{
public:
  virtual bool equal(SignalPtr  &a,SignalPtr  &b) {return false;}
  virtual bool greater(SignalPtr  &a,SignalPtr  &b) {return false;}
  virtual bool less(SignalPtr  &a,SignalPtr  &b) {return true;}
};


class SinglePhaseComparator : public BOOM::Comparator<SignalPtr>
{
protected:
  int phase;
  ContentType contentType;
  DiscreteDistribution &distr;

  virtual double getLengthPenalty(SignalPtr,int position);
  virtual void getScores(SignalPtr,SignalPtr,double &scoreA,double &scoreB);
public:
  SinglePhaseComparator(int phase,ContentType,DiscreteDistribution &);
  virtual bool equal(SignalPtr &,SignalPtr &);
  virtual bool greater(SignalPtr &,SignalPtr &);
  virtual bool less(SignalPtr &,SignalPtr &);
};


class NoncodingComparator : public SinglePhaseComparator
{
public:
  NoncodingComparator(ContentType,DiscreteDistribution&);
  virtual void getScores(SignalPtr,SignalPtr,double &scoreA,double &scoreB);
};



#endif
