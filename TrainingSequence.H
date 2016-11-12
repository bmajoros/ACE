/****************************************************************
 TrainingSequence.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TrainingSequence_H
#define INCL_TrainingSequence_H
#include "BOOM/Sequence.H"


class TrainingSequence : public Sequence
{
public:
  TrainingSequence();
  TrainingSequence(const BOOM::String &seq,Alphabet &);
  int getBoostCount() const;
  void adjustBoostCount(int by);
  virtual Sequence *reverseComplement(Alphabet &);
  virtual void getSubsequence(int begin,int len,Sequence &seq) const;
private:
  int boostCount;
};


#endif
