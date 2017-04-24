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
#include "BOOM/Comparator.H"
#include "BOOM/Interval.H"
using namespace std;
using namespace BOOM;

class LightEdge;

class LightVertex {
public:
  LightVertex(const String &substrate,SignalType,int begin,int end,
	      float score,Strand,int ID);
  virtual ~LightVertex() {}
  virtual void addEdgeIn(LightEdge *);
  virtual void addEdgeOut(LightEdge *);
  virtual SignalType getType() const;
  inline SignalType getSignalType() const { return getType(); }
  inline bool isSpliceSite() { return ::isSpliceSite(type); }
  virtual int getBegin() const;
  virtual int getEnd() const;
  virtual Interval getInterval() const;
  virtual int distanceTo(const LightVertex &) const;
  virtual int getID() const;
  virtual void setID(int);
  virtual float getScore() const;
  virtual void setScore(float);
  virtual Strand getStrand() const;
  virtual Vector<LightEdge*> &getEdgesIn();
  virtual Vector<LightEdge*> &getEdgesOut();
  virtual void printOn(ostream &,String vertexType="vertex");
  virtual bool isSupported() const;
  virtual void setSupport(bool);
  virtual bool isBroken() const { return !supported; }
  virtual bool setBroken(bool b) { supported=!b; }
  virtual void setAnnotated() { annotated=true; }
  virtual void setAnnotated(bool a) { annotated=a; }
  virtual bool isAnnotated() { return annotated; }
  virtual const String &getSubstrate() const { return substrate; }
  virtual void dropEdgeIn(LightEdge *);
  virtual void dropEdgeOut(LightEdge *);
  virtual bool operator==(const LightVertex &) const;
protected:
  Vector<LightEdge*> edgesIn, edgesOut;
  SignalType type;
  float score;
  Strand strand;
  int begin, end;
  int ID;
  String substrate;
  bool supported; // False if broken splice site
  bool annotated; // True if the feature is annotated
};

ostream &operator<<(ostream &,const LightVertex &);

class LightVertexComparator : public Comparator<LightVertex*> {
  bool equal(LightVertex *&a,LightVertex *&b) 
  {return a->getBegin() == b->getBegin();}
  bool greater(LightVertex *&a,LightVertex *&b)
  {return a->getBegin() > b->getBegin();}
  bool less(LightVertex *&a,LightVertex *&b)
  {return a->getBegin() < b->getBegin();}
};

#endif

