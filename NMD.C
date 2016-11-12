/****************************************************************
 NMD.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "NMD.H"
#include "BOOM/CodonIterator.H"
using namespace std;
using namespace BOOM;



ostream &operator<<(ostream &os,NMD_TYPE t)
{
  switch(t) {
  case NMD_NONE:       os<<"none"; break;
  case NMD_NMD:        os<<"NMD"; break;
  case NMD_TRUNCATION: os<<"truncation"; break;
  case NMD_NO_STOP:    os<<"nonstop"; break;
  case NMD_NO_START:   os<<"noncoding"; break;
  }
}


NMD::NMD(int dist)
  : distParm(dist)
{
  // ctor
}



NMD_TYPE NMD::predict(GffTranscript &transcript,const String &substrate,
		      int &distance)
{
  distance=0;
  const int numExons=transcript.getNumExons();
  if(numExons<1) return NMD_NO_START;

  Vector<GffExon*> rawExons; // combines CDS and UTR to get actual exons
  transcript.getRawExons(rawExons);
  const int lastExonLen=rawExons.getLast()->length();

  const int lastEJC=transcript.getSplicedLength()-lastExonLen;
  CodonIterator iter(transcript,rawExons,substrate);
  Codon codon;
  NMD_TYPE ret=NMD_NO_STOP;
  while(iter.nextCodon(codon)) {
    if(codon.isStop()) {
      distance=lastEJC-codon.splicedCoord;
      if(distance>=distParm) ret=NMD_NMD;
      //else if(iter.nextCodon(codon)) ret=NMD_TRUNCATION;
      else ret=NMD_NONE;
      break;
    }
  }
  GffTranscript::deleteExons(rawExons);
  return ret;
}
