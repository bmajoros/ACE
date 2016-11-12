/****************************************************************
 IntronDepthProfile.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_IntronDepthProfile_H
#define INCL_IntronDepthProfile_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "RnaJunctions.H"
using namespace std;
using namespace BOOM;


class IntronDepthProfile {
public:
  IntronDepthProfile(RnaJunctions &,int length);
  float operator[](int) const;
  void rescale(float factor);
  float getMax() const;
protected:
  Array1D<float> A;
  void construct(RnaJunctions &);
  void normalize();
};

#endif

