/****************************************************************
 SignalComparator.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "SignalComparator.H"
#include "Propagator.H"
#include <iostream>
using namespace std;


//==============================================================
//               SinglePhaseComparator methods
//==============================================================

SinglePhaseComparator::SinglePhaseComparator(int phase,
					     ContentType contentType,
					     DiscreteDistribution &distr)
  : phase(phase),
    contentType(contentType),
    distr(distr)
{
}



double SinglePhaseComparator::getLengthPenalty(SignalPtr signal,
					       int position)
{
  int length=position-signal->posOfBaseFollowingConsensus();
  //cout<<"length="<<length<<endl;
  double lengthPenalty=distr.getLogP(length);
  //cout<<" {"<<length<<" "<<signal->posOfBaseFollowingConsensus()<<" "<<position<<"}";
  return lengthPenalty;
}



void SinglePhaseComparator::getScores(SignalPtr a,SignalPtr b,
				      double &scoreA,double &scoreB)
{
  Propagator &propA=a->getPropagator(contentType);
  Propagator &propB=b->getPropagator(contentType);

  int pos=max(propA.getPosition(),propB.getPosition())+1;

  //if(phase<0) throw "BAD";
  scoreA=propA[phase]+getLengthPenalty(a,pos);
  scoreB=propB[phase]+getLengthPenalty(b,pos);
  //cout<<*a<<" "<<*b<<endl;
}



bool SinglePhaseComparator::equal(SignalPtr &a,SignalPtr &b)
{
  double scoreA, scoreB;
  getScores(a,b,scoreA,scoreB);
  return scoreA==scoreB;
}



bool SinglePhaseComparator::greater(SignalPtr &a,SignalPtr &b)
{
  double scoreA, scoreB;
  getScores(a,b,scoreA,scoreB);
  //if(scoreA>scoreB) cout<<"greater "<<scoreA<<" "<<scoreB<<endl; else cout<<"NOT greater "<<scoreA<<" "<<scoreB<<endl;
  return scoreA>scoreB;
}



bool SinglePhaseComparator::less(SignalPtr &a,SignalPtr &b)
{
  //cout<<"less "<<a<<" "<<b<<endl;
  double scoreA, scoreB;
  getScores(a,b,scoreA,scoreB);
  //if(scoreA<scoreB) cout<<"less "<<scoreA<<" "<<scoreB<<endl; else cout<<"NOT less "<<scoreA<<" "<<scoreB<<endl;
  return scoreA<scoreB;
}




//=================================================================
//                  NoncodingComparator methods
//=================================================================

NoncodingComparator::NoncodingComparator(ContentType contentType,
					 DiscreteDistribution &distr)
  : SinglePhaseComparator(-1,contentType,distr)
{
}



void NoncodingComparator::getScores(SignalPtr a,SignalPtr b,
				    double &scoreA,
				    double &scoreB)
{
  Propagator &propA=a->getPropagator(contentType);
  Propagator &propB=b->getPropagator(contentType);

  //cout<<"PROPAGATOR POSITIONS: "<<propA.getPosition()<<" "<<propB.getPosition()<<endl;
  //int pos=max(propA.getPosition(),propB.getPosition())+1;
  if(propA.getPosition()!=propB.getPosition()) throw "UNEQUAL PROPAGATOR POSTIONS IN NoncodingComparator::getScores()";

  int pos=propA.getPosition()+1;
    // both propagators are at the same position!

  int phaseA=(a->getStrand()==FORWARD_STRAND ? 0 : 2);
  int phaseB=(b->getStrand()==FORWARD_STRAND ? 0 : 2);

  scoreA=propA[phaseA]+getLengthPenalty(a,pos);
  scoreB=propB[phaseB]+getLengthPenalty(b,pos);
  //cout<<propA[phaseA]<<":"<<phaseA<<" "<<propB[phaseB]<<":"<<phaseB<<"    "<<propA[2-phaseA]<<"~"<<propB[2-phaseB]<<" "<<*a<<"&"<<*b<<" "<<propA<<"&"<<propB<<" "<<a<<"&"<<b<<"====== "<<getLengthPenalty(a,pos)<<" "<<getLengthPenalty(b,pos)<<endl;
  //cout<<"comparing scores: "<<scoreA<<" "<<scoreB<<" "<<*a<<" "<<*b<<endl;

  //###DEBUGGING:
  //scoreA=propA[phaseA]; scoreB=propB[phaseB];

  //cout<<propA<<" "<<propB<<endl;
}

