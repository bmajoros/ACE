/****************************************************************
 GeometricDistribution.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#include "GeometricDistribution.H"
#include <iostream>
#include <math.h>
#include "BOOM/Exceptions.H"
using namespace std;

GeometricDistribution::GeometricDistribution(int meanLength)
  : usingRatios(false)
{
  //double q=1.0/meanLength;
  double q=1.0/(meanLength+1); // ### 4/19/2011
  logQ=log(q);
  logOneMinusQ=log(1-q);
  logMean=getLogP(meanLength);
}



double GeometricDistribution::getLogP(unsigned len)
{
  // log[ (1-q)^(len-1) * q ] == (len-1)*log(1-q)+log(q)

  if(len<0) INTERNAL_ERROR;
  /*
  if(len<1) return logQ;
  return (len-1)*logOneMinusQ+logQ;//-(usingRatios ? logMean : 0);
  */
  return len*logOneMinusQ+logQ;// ### 4/19/2011
}



void GeometricDistribution::useLogLikelihoodRatios()
{
  usingRatios=true;
}




