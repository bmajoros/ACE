/****************************************************************
 Signal.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Signal_H
#define INCL_Signal_H
#include <iostream>
#include "BOOM/List.H"
#include "BOOM/Comparator.H"
#include "BOOM/Array1D.H"
#include "genezilla.H"
#include "SignalType.H"
class Edge;
class SignalSensor;
class Propagator;
class GarbageCollector;

class Signal : public PointedTo
{
#ifdef EXPLICIT_GRAPHS
  BOOM::Set<Edge*> edgesIn;
  BOOM::Set<Edge*> edgesOut;
  int index; // index into graph's vertex array
#endif

  BOOM::SmartPointer<Signal> predecessors[3]; // optimal predecessors, one in each frame
  int contextWindowPosition;
  Array1D<Propagator*> propagators;
  SignalType signalType;
  SignalSensor &sensor;
  double signalScore; //score for the entire context window (not just cons.)
  double inductiveScores[3]; // includes signalScore!
  bool annotated; // is annotated as a true site

  void initializePropagators(double signalScore);
public:
  Signal(int contextWindowPosition,double signalScore,SignalSensor &,
	 GarbageCollector &,SignalType);
  virtual ~Signal();
  void setAnnotated() { annotated=true; }
  bool isAnnotated() { return annotated; }
  void dropSignalScores();
  void dropContentScores();
  Propagator &getPropagator(ContentType);
  void setPredecessor(int phase,BOOM::SmartPointer<Signal>);
  BOOM::SmartPointer<Signal> getPredecessor(int phase);
  int getContextWindowPosition();
  int getContextWindowEnd(); // index of first base after the window
  int getConsensusPosition() const;
  virtual int getContextWindowLength();
  virtual int getConsensusLength();
  int frameOfBaseFollowingConsensus();
  int posOfBaseFollowingConsensus();
  BOOM::Set<ContentType> &belongsInWhichQueues();
  BOOM::Set<ContentType> &linksBackToWhichQueues();
  Strand getStrand();
  SignalType getSignalType() const;
  double contextWindowScore(); // signal score, not the inductive score
  SignalSensor &getSensor();
  double &getInductiveScore(int phase);
  double priorInductiveScore(int phase); // not including signal score
  double posteriorInductiveScore(int phase); // including signal score
  double precedingFeatureScore(int phase); // =content+length+trans scores
  bool isLeftTerminus();
  void setSignalType(SignalType);
  void printOn(ostream &) const;

#ifdef EXPLICIT_GRAPHS
  bool isRightTerminus();
  void addEdgeIn(Edge *);
  void addEdgeOut(Edge *);
  BOOM::Set<Edge*> &getEdgesIn();
  BOOM::Set<Edge*> &getEdgesOut();
  void dropEdgeIn(Edge*);
  void dropEdgeOut(Edge*);
  inline int getIndex() {return index;}
  inline void setIndex(int x) {index=x;}
  Edge *findEdgeInFrom(BOOM::SmartPointer<Signal> pred);
#endif
};

ostream &operator<<(ostream &,const Signal &);


struct SignalPosComparator : public Comparator<Signal*> {
  bool equal(Signal* &a,Signal* &b);
  bool greater(Signal* &a,Signal* &b);
  bool less(Signal* &a,Signal* &b);
};





#endif
