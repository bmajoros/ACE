/****************************************************************
 EmpiricalDistribution.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_EmpiricalDistribution_H
#define INCL_EmpiricalDistribution_H

#include "BOOM/String.H"
#include "BOOM/Vector.H"
#include "DiscreteDistribution.H"

class EmpiricalDistribution : public DiscreteDistribution
{
  typedef pair<unsigned,double> EmpiricalDistributionElement;
  BOOM::Vector<EmpiricalDistributionElement*> v;
  double smallestElemLogP, largestElemLogP, mean;
  int binSize;
  bool useInterpolation;

  unsigned binarySearch(unsigned elem);
  void load(const BOOM::String &filename);
  void normalize();
  inline double interpolate(unsigned x1,double y1,unsigned x2,double y2,
			    unsigned x);
public:
  EmpiricalDistribution(const BOOM::String &filename,
			bool shouldInterpolate=false);
  virtual ~EmpiricalDistribution();
  virtual double getLogP(unsigned);
  virtual void useLogLikelihoodRatios();
};


#endif
