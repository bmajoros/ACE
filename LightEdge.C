/****************************************************************
 LightEdge.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "LightEdge.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
#include "LightVertex.H"
using namespace std;
using namespace BOOM;

LightEdge::LightEdge(const String &substrate,ContentType t,LightVertex *l,
		     LightVertex *r,int begin,int end,Strand s,int ID)
  : type(t), left(l), right(r), begin(begin), end(end), strand(s),
    substrate(substrate), ID(ID), supported(false)
{
  score[0]=score[1]=score[2]=NEGATIVE_INFINITY;
}



ContentType LightEdge::getType() const
{
  return type;
}



int LightEdge::getBegin() const
{
  return begin;
}



int LightEdge::getEnd() const
{
  return end;
}



LightVertex *LightEdge::getLeft()
{
  return left;
}



LightVertex *LightEdge::getRight()
{
  return right;
}



float LightEdge::getScore(int frame) const
{
  return score[frame];
}



void LightEdge::setScore(int frame,float s)
{
  score[frame]=s;
}



Strand LightEdge::getStrand() const
{
  return strand;
}



void LightEdge::printOn(ostream &os)
{
  int supportString=supported ? 1 : 0;
  bool useFrames=::isIntron(type) || ::isCoding(type);
  //bool needNewline=false;
  if(useFrames) {
    for(int i=0 ; i<3 ; ++i)
      if(isFinite(score[i])) {
	//if(needNewline) os<<endl;
	os<<substrate<<"\tedge\t"<<contentTypeNiceString(type)<<"\t"<<begin+1
	  <<"\t"<<end<<"\t"<<score[i]<<"\t"<<strand<<"\t"<<i
	  <<"\tleft="<<left->getID()<<"; right="<<right->getID()
	  <<"; edgeID="<<ID<<"; sup="<<supportString<<";"
	  <<endl;
	//needNewline=true;
      }
  }
  else
    os<<substrate<<"\tedge\t"<<contentTypeNiceString(type)<<"\t"<<begin+1
      <<"\t"<<end<<"\t"<<score[0]<<"\t"<<strand<<"\t.\tleft="
      <<left->getID()<<"; right="<<right->getID()<<"; edgeID="<<ID<<";"
      <<" sup="<<supportString<<";"
      <<endl;
}



ostream &operator<<(ostream &os,const LightEdge &e)
{
  e.printOn(os);
  return os;
}



Array1D<int> LightEdge::getFrames() const
{
  int n=0;
  for(int i=0 ; i<3 ; ++i) if(isFinite(score[i])) ++n;
  Array1D<int> A(n);
  for(int i=0, j=0 ; i<3 ; ++i) if(isFinite(score[i])) A[j++]=i;
  return A;
}


void LightEdge::subsumeVertexScores()
{
  float leftScore=left->getScore();
  float rightScore=right->getScore();
  float leftPlusRight=leftScore+rightScore;
  for(int i=0 ; i<3 ; ++i) 
    if(isFinite(score[i]))
      score[i]+=leftPlusRight;
}



bool LightEdge::isCoding() const
{
  return ::isCoding(type);
}



bool LightEdge::isExon() const
{
  return ::isUTR(type) || ::isCoding(type);
}



bool LightEdge::isIntron() const
{
  return ::isIntron(type);
}



bool LightEdge::isIntergenic() const
{
  return ::isIntergenic(type);
}



int LightEdge::getLength() const
{
  return end-begin;
}



int LightEdge::propagateForward(int phase) const
{
  if(isIntergenic()) return (right->getStrand()==FORWARD_STRAND ? 0 : 2);
  if(!isCoding()) return phase;
  int length=getEnd()-getBegin();
  switch(getStrand())
    {
    case FORWARD_STRAND:
      return (phase+length)%3;
    case REVERSE_STRAND:
      return posmod(phase-length);
    }
  INTERNAL_ERROR;
}



int LightEdge::propagateBackward(int phase) const
{
  if(isIntergenic()) return (left->getStrand()==FORWARD_STRAND ? 0 : 2);
  if(!isCoding()) return phase;
  int length=getEnd()-getBegin();
  switch(getStrand())
    {
    case FORWARD_STRAND:
      return posmod(phase-length);
    case REVERSE_STRAND:
      return (phase+length)%3;
    }
}



int LightEdge::posmod(int x) const
{
  int f=x%3;
  if(f>=0) return f;
  return f+3;
}



bool LightEdge::isSupported() const
{
  return supported;
}



void LightEdge::setSupport(bool s)
{
  supported=s;
}



