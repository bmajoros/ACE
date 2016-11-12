/****************************************************************
 StartCodonFinder.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "StartCodonFinder.H"
#include "BOOM/Sequence.H"
using namespace std;
using namespace BOOM;



int StartCodonFinder::findStartCodon(const GffTranscript &transcript,
				     const Vector<GffExon*> &rawExons,
				     const String &genome,int pos,
				     const SignalSensors &sensors)
{
  transcript.loadSequence(genome);
  String splicedTranscript=transcript.getFullSequence();
  Sequence splicedSeq(splicedTranscript,alphabet);
  const int splicedLen=splicedTranscript.length();
  pos=GffTranscript::genomicToSplicedCoords(pos,rawExons);
  SignalSensor *sensor=sensors.startCodonSensor;
  SignalSensor *shortSensor=sensors.shortStartSensor;
  const float cutoff=sensor->getCutoff();
  const int offset=sensor->getConsensusOffset();
  const int windowLen=sensor->getContextWindowLength();
  bool firstPos=true;
  while(pos-offset+windowLen<splicedLen) {
    const float score=
      pos-offset<0 || firstPos ?
      shortSensor->getLogP(splicedSeq,splicedTranscript,pos) :
      sensor->getLogP(splicedSeq,splicedTranscript,pos-offset);
    if(score>=cutoff)
      return GffTranscript::splicedToGenomicCoords(pos,rawExons);
    ++pos; firstPos=false;
  }
  return -1;
}
