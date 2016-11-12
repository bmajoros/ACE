/****************************************************************
 TataCapSignal.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TataCapSignal_H
#define INCL_TataCapSignal_H
#include "Signal.H"

class TataCapSignal : public Signal
{
public:
  TataCapSignal(GarbageCollector &,int contextWindowPosition,
		int contextWindowLength,
		int consensusLength,double signalScore,SignalSensor &);
  virtual int getContextWindowLength();
  virtual int getConsensusLength();
private:
  int contextWindowLength;
  int consensusLength; // bogus -- but doesn't matter (too much)
};


#endif
