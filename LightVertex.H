/****************************************************************
 LightVertex.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_LightVertex_H
#define INCL_LightVertex_H
#include <iostream>
#include "BOOM/Vector.H"
#include "SignalType.H"
#include "BOOM/Strand.H"
using namespace std;
using namespace BOOM;

class LightEdge;

class LightVertex {
public:
  LightVertex(const String &substrate,SignalType,int begin,int end,
	      float score,Strand,int ID);
  void addEdgeIn(LightEdge *);
  void addEdgeOut(LightEdge *);
  SignalType getType() const;
  inline SignalType getSignalType() const { return getType(); }
  int getBegin() const;
  int getEnd() const;
  int getID() const;
  float getScore() const;
  void setScore(float);
  Strand getStrand() const;
  Vector<LightEdge*> &getEdgesIn();
  Vector<LightEdge*> &getEdgesOut();
  void printOn(ostream &,String vertexType="vertex");
  bool isSupported() const;
  void setSupport(bool);
  void setAnnotated() { annotated=true; }
  bool isAnnotated() { return annotated; }
protected:
  Vector<LightEdge*> edgesIn, edgesOut;
  SignalType type;
  float score;
  Strand strand;
  int begin, end;
  int ID;
  String substrate;
  bool supported, annotated;
};

ostream &operator<<(ostream &,const LightVertex &);

#endif

