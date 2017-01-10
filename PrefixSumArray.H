/****************************************************************
 PrefixSumArray.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_PrefixSumArray_H
#define INCL_PrefixSumArray_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "BOOM/Sequence.H"
#include "BOOM/String.H"
#include "ContentSensor.H"
using namespace std;
using namespace BOOM;


class PrefixSumArray {
public:
  PrefixSumArray(int length=0);
  void resize(int length);
  void compute(const ContentSensor &,const Sequence &,const String&);
  double getInterval(int begin,int end) const;
  double operator[](int) const;
protected:
  Array1D<double> A;
};

#endif

