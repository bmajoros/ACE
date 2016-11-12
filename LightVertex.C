/****************************************************************
 LightVertex.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "LightVertex.H"
using namespace std;
using namespace BOOM;

LightVertex::LightVertex(const String &substrate,SignalType t,int begin,
			 int end,float score,Strand s,int ID)
  : type(t), begin(begin), end(end), score(score), strand(s), ID(ID),
    substrate(substrate), supported(false)
{
  // ctor
}



void LightVertex::addEdgeIn(LightEdge *e)
{
  edgesIn.push_back(e);
}



void LightVertex::addEdgeOut(LightEdge *e)
{
  edgesOut.push_back(e);
}



SignalType LightVertex::getType() const
{
  return type;
}



int LightVertex::getBegin() const
{
  return begin;
}



int LightVertex::getEnd() const
{
  return end;
}



float LightVertex::getScore() const
{
  return score;
}



void LightVertex::setScore(float s)
{
  score=s;
}



Strand LightVertex::getStrand() const
{
  return strand;
}



Vector<LightEdge*> &LightVertex::getEdgesIn()
{
  return edgesIn;
}



Vector<LightEdge*> &LightVertex::getEdgesOut()
{
  return edgesOut;
}



void LightVertex::printOn(ostream &os,String vertexType)
{
  int supportString=supported ? 1 : 0;
  os<<substrate<<"\t"<<vertexType<<"\t"<<signalTypeToName(type)
    <<"\t"<<begin+1<<"\t"<<
    end<<"\t"<<score<<"\t"<<strand<<"\t.\tID="<<ID<<"; in="<<edgesIn.size()
    <<"; out="<<edgesOut.size()<<"; sup="<<supportString<<";";
}



int LightVertex::getID() const
{
  return ID;
}



ostream &operator<<(ostream &os,const LightVertex &v)
{
  v.printOn(os);
  return os;
}



bool LightVertex::isSupported() const
{
  return supported;
}



void LightVertex::setSupport(bool s)
{
  supported=s;
}



