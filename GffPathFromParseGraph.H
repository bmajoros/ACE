/****************************************************************
 GffPathFromParseGraph.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_GffPathFromParseGraph_H
#define INCL_GffPathFromParseGraph_H

#include "BOOM/Vector.H"
#include "BOOM/GffTranscript.H"
#include "BOOM/String.H"
#include "BOOM/Sequence.H"
#include "Signal.H"
#include "genezilla.H"
#ifdef EXPLICIT_GRAPHS

// This is required by the standards committee screwed up vector<bool>...
typedef char BOOL;

class GeneZilla;
class GffPathFromParseGraph
{
  GeneZilla &genezilla;

  void exonCoordsToSigPos(ContentType,int begin,int end,int &pos1,
			  int &pos2);
  SignalPtr getSignal(SignalType signalType,int pos,const Sequence &seq,
		      const BOOM::String &seqStr,BOOL &found);
public:
  GffPathFromParseGraph(GeneZilla &);
  BOOM::Vector<SignalPtr> *getPathFromGff(
				   BOOM::Vector<BOOM::GffTranscript*> &,
					  const Sequence &,
					  const BOOM::String &seqStr,
					  BOOM::Vector<BOOL> &found);
  BOOM::Vector<int> *getSignalCoordinates(
				   BOOM::Vector<BOOM::GffTranscript*> &);
};


#endif
#endif
