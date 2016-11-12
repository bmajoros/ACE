/****************************************************************
 TataCapSignal.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "TataCapSignal.H"
#include <iostream>


TataCapSignal::TataCapSignal(GarbageCollector &gc,
			     int contextWindowPosition,
			     int contextWindowLength,
			     int consensusLength,double signalScore,
			     SignalSensor &signalSensor)
  : Signal(contextWindowPosition,signalScore,signalSensor,gc),
    consensusLength(consensusLength),
    contextWindowLength(contextWindowLength)
{
}



int TataCapSignal::getContextWindowLength()
{
  return contextWindowLength;
}



int TataCapSignal::getConsensusLength()
{
  return consensusLength;
}


