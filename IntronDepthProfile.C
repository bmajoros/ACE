/****************************************************************
 IntronDepthProfile.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "IntronDepthProfile.H"
using namespace std;
using namespace BOOM;

IntronDepthProfile::IntronDepthProfile(RnaJunctions &junctions,int length)
  : A(length)
{
  construct(junctions);
  //normalize();
}



float IntronDepthProfile::operator[](int i) const
{
  return A[i];
}



void IntronDepthProfile::construct(RnaJunctions &junctions)
{
  int L=A.size();
  A.setAllTo(0);
  int numJ=junctions.getNumJunctions();
  for(int i=0 ; i<numJ ; ++i) {
    const RnaJunction &junction=junctions[i];
    int begin=junction.getBegin(), end=junction.getEnd();
    float depth=junction.getDepth();
    for(int j=begin ; j<end ; ++j) A[j]+=depth;
  }
}



void IntronDepthProfile::normalize()
{
  int L=A.size();
  float max=0;
  for(int i=0 ; i<L ; ++i) {
    float depth=A[i];
    if(depth>max) max=depth;
  }
  for(int i=0 ; i<L ; ++i) A[i]/=max;
}



void IntronDepthProfile::rescale(float factor)
{
  int L=A.size();
  for(int i=0 ; i<L ; ++i)
    A[i]*=factor;
}



float IntronDepthProfile::getMax() const
{
  float m=0;
  int L=A.size();
  for(int i=0 ; i<L ; ++i) {
    float x=A[i];
    if(x>m) m=x;
  }
  return m;
}

