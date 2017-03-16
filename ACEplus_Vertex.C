/****************************************************************
 ACEplus_Vertex.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ACEplus_Vertex.H"
using namespace std;
using namespace BOOM;

ACEplus_Vertex::ACEplus_Vertex(const String &substrate,SignalType t,int begin,
			       int end,double score,Strand s,int ID)
  : LightVertex(substrate,t,begin,end,score,s,ID),
    threshold(0.0), refScore(0.0)
{
  // ctor
}



void ACEplus_Vertex::setRefScore(double s)
{
  refScore=s;
}



void ACEplus_Vertex::setThreshold(double t)
{
  threshold=t;
}



void ACEplus_Vertex::setSeq(const String &s)
{
  seq=s;
}



double ACEplus_Vertex::getRefScore() const
{
  return refScore;
}



double ACEplus_Vertex::getThreshold() const
{
  return threshold;
}



const String &ACEplus_Vertex::getSeq() const
{
  return seq;
}



void ACEplus_Vertex::setRawScore(double s)
{
  rawScore=s;
}



double ACEplus_Vertex::getRawScore() const
{
  return rawScore;
}




