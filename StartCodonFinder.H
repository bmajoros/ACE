/****************************************************************
 StartCodonFinder.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_StartCodonFinder_H
#define INCL_StartCodonFinder_H
#include <iostream>
#include "BOOM/GffTranscript.H"
#include "BOOM/Vector.H"
#include "BOOM/String.H"
#include "SignalSensors.H"

using namespace std;
using namespace BOOM;


class StartCodonFinder {
public:
  static int findStartCodon(const GffTranscript &transcript,
			    const Vector<GffExon*> &rawExons,
			    const String &genome,int pos,
			    const SignalSensors &);
};

#endif

