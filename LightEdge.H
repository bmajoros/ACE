/****************************************************************
 LightEdge.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_LightEdge_H
#define INCL_LightEdge_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Strand.H"
#include "ContentType.H"
using namespace std;
using namespace BOOM;

class LightVertex;

class LightEdge {
public:
  LightEdge(const String &substrate,ContentType,LightVertex *,LightVertex *,
	    int begin,int end,Strand,int ID);
  ContentType getType() const;
  inline ContentType getEdgeType() const { return getType(); }
  bool isExon() const;
  bool isCoding() const;
  bool isIntron() const;
  bool isIntergenic() const;
  int getBegin() const;
  int getEnd() const;
  int getLength() const;
  LightVertex *getLeft();
  LightVertex *getRight();
  float getScore(int frame) const;
  void setScore(int frame,float score);
  Strand getStrand() const;
  void printOn(ostream &);
  Array1D<int> getFrames() const;
  void subsumeVertexScores();
  int propagateForward(int phase) const;
  int propagateBackward(int phase) const;
  bool isSupported() const;
  void setSupport(bool);
  int getID() const { return ID; }
protected:
  LightVertex *left, *right;
  ContentType type;
  int begin, end;
  float score[3];
  int ID;
  Strand strand;
  String substrate;
  bool supported;
  int posmod(int) const;
};

ostream &operator<<(ostream &,const LightEdge &);

struct EdgeEndComparator : public Comparator<LightEdge*> {
  bool equal(LightEdge* &a,LightEdge* &b)   {return a->getEnd()==b->getEnd();}
  bool greater(LightEdge* &a,LightEdge* &b) {return a->getEnd()>b->getEnd();}
  bool less(LightEdge* &a,LightEdge* &b)    {return a->getEnd()<b->getEnd();}
};


#endif

