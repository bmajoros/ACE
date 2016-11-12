/****************************************************************
 TrainingSequence.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "TrainingSequence.H"
#include <iostream>


TrainingSequence::TrainingSequence(const BOOM::String &seq,
				   Alphabet &alphabet)
  : boostCount(1),
    Sequence(seq,alphabet)
{
}



TrainingSequence::TrainingSequence()
  : boostCount(1)
{
}



int TrainingSequence::getBoostCount() const
{
  return boostCount;
}



void TrainingSequence::adjustBoostCount(int by)
{
  boostCount+=by;
}



Sequence *TrainingSequence::reverseComplement(Alphabet &alphabet)
{
  TrainingSequence *rc=new TrainingSequence;
  Sequence::reverseComplement(alphabet,*rc);
  rc->boostCount=boostCount;
  return rc;
}



void TrainingSequence::getSubsequence(int begin,int len,Sequence &seq) const
{
  Sequence::getSubsequence(begin,len,seq);
  TrainingSequence *ts=dynamic_cast<TrainingSequence*>(&seq);
  if(ts) ts->boostCount=boostCount;
}
