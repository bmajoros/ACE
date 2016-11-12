/****************************************************************
 DiscreteDistribution.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_DiscreteDistribution_H
#define INCL_DiscreteDistribution_H

enum DistributionType
  {
    EMPIRICAL_DISTRIBUTION,
    GEOMETRIC_DISTRIBUTION
  };

class DiscreteDistribution
{
public:
  DiscreteDistribution() {}
  virtual ~DiscreteDistribution() {}
  virtual double getLogP(unsigned)=0;
  virtual void useLogLikelihoodRatios()=0;
};


#endif
