/****************************************************************
 SignalPrinter.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SignalPrinter.H"
using namespace std;
using namespace BOOM;


String SignalPrinter::print(SignalSensor &sensor,
			    int windowPos,const String &genome)
{
  String window;
  const int offset=sensor.getConsensusOffset();
  const int windowLen=sensor.getContextWindowLength();
  const int consensusLen=sensor.getConsensusLength();
  window=
    genome.substring(windowPos,offset).tolower()
    + "_"
    + genome.substring(windowPos+offset,consensusLen)
    + "_"
    + genome.substring(windowPos+offset+consensusLen,
		       windowLen-offset-consensusLen).tolower();
  return window;
}


