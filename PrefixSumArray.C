/****************************************************************
 PrefixSumArray.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "PrefixSumArray.H"
using namespace std;
using namespace BOOM;

PrefixSumArray::PrefixSumArray(int length)
  : A(length)
{
  // ctor
}


void PrefixSumArray::resize(int length)
{
  A.resize(length);
}



void PrefixSumArray::compute(const ContentSensor &sensor,const Sequence &seq,
			     const String &str)
{
  const int L=A.size();
  if(L>seq.getLength() || L>str.getLength()) INTERNAL_ERROR;
  double sum=0;
  for(int i=0 ; i<L ; ++i) {
    sum+=sensor.scoreSingleBase(seq,str,i,seq[i],str[i]);
    A[i]=sum;
  }
}



double PrefixSumArray::getInterval(int begin,int end) const
{
  if(begin>=end) return 0.0;
  double beginScore=begin>0 ? A[begin-1] : 0.0;
  return A[end-1]-beginScore;
}



double PrefixSumArray::operator[](int i) const
{
  return A[i];
}


