/****************************************************************
 TrellisLink.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "TrellisLink.H"
using namespace std;
using namespace BOOM;

TrellisLink::TrellisLink(TrellisLink *pred,LightEdge *edge,int phase,
			 float score)
  : pred(pred), edge(edge), phase(phase), score(score)
{
}



int TrellisLink::getPhase() const
{
  return phase;
}



TrellisLink *TrellisLink::getPred() const
{
  return pred;
}



LightEdge *TrellisLink::getEdge() const
{
  return edge;
}



float TrellisLink::getScore() const
{
  return score;
}



void TrellisLink::setScore(float s)
{
  score=s;
}



