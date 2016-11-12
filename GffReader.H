/****************************************************************
 GffReader.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_GffReader_H
#define INCL_GffReader_H
#include "BOOM/String.H"
#include "BOOM/Map.H"
#include "BOOM/Vector.H"
#include "SignalSensor.H"
#include "GarbageCollector.H"
#include "ParseGraph.H"
#include "genezilla.H"
class GeneZilla;

class GffReader
{
public:
  typedef BOOM::Vector<SignalPtr> SignalList;
  GffReader(const BOOM::String &filename,BOOM::Vector<SignalSensor*> &,
	    GarbageCollector &);
  ~GffReader();
#ifdef EXPLICIT_GRAPHS
  ParseGraph *toGraph(GeneZilla &);
#endif
  SignalList &getSignals();
  SignalList *signalsAt(int position);
  SignalPtr getSignal(int position,SignalType);
  SignalPtr getSignal(int position,SignalType,SignalPtr &pred);
  void addLeftTerminus(SignalPtr );
  void addRightTerminus(SignalPtr );
private:
  SignalList signals;
  BOOM::Map<int,SignalList*> posToSignals;
  BOOM::Map<SignalType,SignalSensor*> signalSensors;
  GarbageCollector &gc;

  void load(const BOOM::String &);
  void newSignal(int contextWindowPos,float score,SignalSensor &);
  void linkBack(SignalPtr);
  void createLeftTerminus(SignalType);
  int getPhase(SignalPtr pred);
};


#endif
