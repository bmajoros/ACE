/**************************************************************
 EmpiricalDistribution.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include "EmpiricalDistribution.H"
#include "BOOM/Random.H"
#include <math.h>
#include <iostream>
#include <fstream>


EmpiricalDistribution::EmpiricalDistribution(const BOOM::String &filename,
					     bool shouldInterpolate)
  : useInterpolation(shouldInterpolate)
{
  // ctor

  load(filename);
  normalize();

  smallestElemLogP=log(v[0]->second);
  largestElemLogP=log(v[v.size()-1]->second);
}



EmpiricalDistribution::~EmpiricalDistribution()
{
  BOOM::Vector<EmpiricalDistributionElement*>::iterator cur=v.begin(),
    end=v.end();
  for(; cur!=end ; ++cur) delete *cur;
}



void EmpiricalDistribution::normalize()
{
  // Sum the values (even those not explicitly represented)
  int n=v.size(), halfBinSize=binSize/2;
  double sum=0;
  for(int i=0 ; i<n ; ++i)
    {
      EmpiricalDistributionElement &elem=*v[i];
      int x=elem.first+halfBinSize;
      double y=elem.second;
      sum+=y*binSize; // *binSize is necessary so all probabilities are
                      // comparable, even across different histograms and
                      // feature types (i.e., coding & noncoding)
    }

  // Divide by the sum
  for(int i=0 ; i<n ; ++i) v[i]->second/=sum;
}



unsigned EmpiricalDistribution::binarySearch(unsigned elem)
{
  unsigned begin=0, end=v.size();
  while(begin<end)
    {
      unsigned mid=unsigned((begin+end)/2);
      unsigned midElem=v[mid]->first;
      if(elem>midElem) begin=mid+1;
      else end=mid;
    }
  return begin;
}



void EmpiricalDistribution::load(const BOOM::String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file ")+filename+
		   " in EmpiricalDistribution::load()";
  while(!is.eof())
    {
      unsigned x;
      double y;
      is >> x;
      if(is.eof()) break;
      is >> y;
      v.push_back(new EmpiricalDistributionElement(x,y));
    }
  binSize=v[1]->first-v[0]->first;
}



double EmpiricalDistribution::getLogP(unsigned x)
{
  // Perform binary search to find nearest neighbor
  unsigned index=binarySearch(x);
  if(index<0) return smallestElemLogP;
  if(index>=v.size()) return largestElemLogP;
  EmpiricalDistributionElement &elem=*v[index];
  unsigned foundX=elem.first, x1, x2;
  double foundY=elem.second, y1, y2;
  if(x==foundX) return log(foundY); // ### this log could be cached

  if(useInterpolation)
    {
      // Perform linear interpolation between two surrounding neighbors
      if(x<foundX)
	{
	  if(index==0) return log(interpolate(0,0,foundX,foundY,x));
	  EmpiricalDistributionElement &prevElem=*v[index-1];
	  x1=prevElem.first;
	  y1=prevElem.second;
	  x2=foundX;
	  y2=foundY;
	}
      else // x>foundX
	{
	  if(index==v.size()-1) return largestElemLogP;
	  EmpiricalDistributionElement &nextElem=*v[index+1];
	  x1=foundX;
	  y1=foundY;
	  x2=nextElem.first;
	  y2=nextElem.second;
	}
      double y=interpolate(x1,y1,x2,y2,x);
      return log(y);
    }
  else return log(foundY);
}



double EmpiricalDistribution::interpolate(unsigned x1,double y1,unsigned x2,
				 double y2,unsigned x)
{
  return (y2-y1)/double(x2-x1)*(x-x1)+y1;
}



void EmpiricalDistribution::useLogLikelihoodRatios()
{
  int n=v.size();
  double mean=0;
  double debug=0;
  for(int i=0 ; i<n ; ++i) 
    {
      int x=v[i]->first;
      double y=v[i]->second;
      debug+=y*binSize;
      mean+=(x+binSize/2)*y*binSize;
    }
  double Pmean=exp(getLogP(int(mean)));

  for(int i=0 ; i<n ; ++i) 
    {
      v[i]->second/=Pmean;
    }
}


