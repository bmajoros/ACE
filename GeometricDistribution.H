/****************************************************************
 GeometricDistribution.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_GeometricDistribution_H
#define INCL_GeometricDistribution_H

#include "DiscreteDistribution.H"

class GeometricDistribution : public DiscreteDistribution
{
  double logQ, logOneMinusQ, logMean;
  bool usingRatios;
public:
  GeometricDistribution(int meanLength);
  virtual double getLogP(unsigned);
  virtual void useLogLikelihoodRatios();
};


#endif
