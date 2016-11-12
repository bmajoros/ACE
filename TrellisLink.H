/****************************************************************
 TrellisLink.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TrellisLink_H
#define INCL_TrellisLink_H
#include <iostream>
#include "BOOM/Comparator.H"
#include "LightEdge.H"
using namespace std;
using namespace BOOM;


class TrellisLink {
public:
  TrellisLink(TrellisLink * =NULL,LightEdge * =NULL,int phase=0,
	      float score=0.0);
  TrellisLink *getPred() const;
  LightEdge *getEdge() const;
  float getScore() const;
  void setScore(float);
  int getPhase() const;
protected:
  TrellisLink *pred;
  LightEdge *edge;
  float score;
  int phase;
};


struct TrellisLinkComparator : public Comparator<TrellisLink> {
  virtual bool equal(TrellisLink &a,TrellisLink &b) 
    {return a.getScore()==b.getScore();}
  virtual bool greater(TrellisLink &a,TrellisLink &b) 
    {return a.getScore()>b.getScore();}
  virtual bool less(TrellisLink &a,TrellisLink &b) 
    {return a.getScore()<b.getScore();}
};


#endif

