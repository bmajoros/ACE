/****************************************************************
 GffPathFromParseGraph.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "GffPathFromParseGraph.H"
#include "genezilla.H"
#ifdef EXPLICIT_GRAPHS
#include "GZilla.H"
#include "ParseGraph.H"
#include <iostream>
using namespace std;
#include "BOOM/VectorSorter.H"


GffPathFromParseGraph::GffPathFromParseGraph(GeneZilla &genezilla)
  : genezilla(genezilla)
{
  // ctor
}



SignalPtr GffPathFromParseGraph::getSignal(SignalType signalType,
					   int pos,const Sequence &seq,
					   const BOOM::String &seqStr,
					   BOOL &found)
{
  // First, see if the signal is in the parse graph (and therefore 
  // already scored)
  ParseGraph &graph=genezilla.getParseGraph();
  int index=graph.findSignal(signalType,pos);
  if(index>-1) 
    {
      found=true;
      return graph.getIthVertex(index);
    }

  // Otherwise, we have to invoke the appropriate signal sensor
  // to score this putative signal
  throw "GffPathFromParseGraph::getSignal() needs to be updated";
  /*
  SignalSensor &sensor=genezilla.getSignalSensor(signalType);
  double logP=sensor.getLogP(seq,seqStr,pos);
  int contextWindowPos=pos-sensor.getConsensusOffset();
  found=false;
  SignalPtr p=new Signal(contextWindowPos,logP,sensor,genezilla.getGC());
  p->setIndex(-1);
  return p;
  */
}



BOOM::Vector<SignalPtr> *GffPathFromParseGraph::getPathFromGff(
  BOOM::Vector<BOOM::GffTranscript*> &transcripts,const Sequence &seq,
  const BOOM::String &seqStr,BOOM::Vector<BOOL> &found)
{
  BOOM::Vector<SignalPtr> &path=*new BOOM::Vector<SignalPtr>;
  int numTranscripts=transcripts.size();
  for(int i=0 ; i<numTranscripts ; ++i) {
    BOOM::GffTranscript &transcript=*transcripts[i];
    if(transcript.getSubstrate()!=genezilla.getSubstrateId()) continue;
    int numExons=transcript.getNumExons();
    for(int j=0 ; j<numExons ; ++j) {
      BOOM::GffExon &exon=transcript.getIthExon(j);
      Strand strand=exon.getStrand();
      ContentType contentType=
	::exonTypeToContentType(exon.getExonType(),strand);
      Set<SignalType> leftSignals=::leftSignals(contentType);
      Set<SignalType> rightSignals=::rightSignals(contentType);
      if(leftSignals.size()!=1 || rightSignals.size()!=1) INTERNAL_ERROR;
      SignalType leftSignal=leftSignals.getSingleElement();
      SignalType rightSignal=rightSignals.getSingleElement();
      int begin=exon.getBegin(), end=exon.getEnd(), pos1, pos2;
      exonCoordsToSigPos(contentType,begin,end,pos1,pos2);
      found.push_back(false);
      SignalPtr signal1=
	getSignal(leftSignal,pos1,seq,seqStr,found[found.size()-1]);
      path.push_back(signal1);
      found.push_back(false);
      SignalPtr signal2=
	getSignal(rightSignal,pos2,seq,seqStr,found[found.size()-1]);
      path.push_back(signal2);
    }
  }

  class SignalCmp : public BOOM::Comparator<SignalPtr>
  {
    bool equal(SignalPtr &a,SignalPtr &b) 
      {return a->getConsensusPosition()==b->getConsensusPosition();}
    bool less(SignalPtr &a,SignalPtr &b) 
      {return a->getConsensusPosition()<b->getConsensusPosition();}
    bool greater(SignalPtr &a,SignalPtr &b) 
      {return a->getConsensusPosition()>b->getConsensusPosition();}
  } cmp;
  DirectComparator<BOOL> cmp2;
  BOOM::VectorSorter<SignalPtr> sorter(path,cmp);
  BOOM::VectorSorter<BOOL> boolSorter(found,cmp2);
  BOOM::Vector<int> *indices=sorter.sortAscendByIndex();
  sorter.applyIndexMap(*indices,path);
  boolSorter.applyIndexMap(*indices,found);
  delete indices;
  return &path;
}



void GffPathFromParseGraph::exonCoordsToSigPos(ContentType contentType,
					       int begin,int end,
					       int &pos1,int &pos2)
{
  switch(contentType)
    {
    case INITIAL_EXON:             // ATG - GT
      pos1=begin; pos2=end; break;
    case INTERNAL_EXON:            // AG - GT
      pos1=begin-2; pos2=end; break;
    case FINAL_EXON:               // AG - TAG
      pos1=begin-2; pos2=end-3; break;
    case SINGLE_EXON:              // ATG - TAG
      pos1=begin; pos2=end-3; break;
    case NEG_INITIAL_EXON:         // NEG_GT - NEG_ATG
      pos1=begin-2; pos2=end-3; break;
    case NEG_INTERNAL_EXON:        // NEG_GT - NEG_AG
      pos1=begin-2; pos2=end; break;
    case NEG_FINAL_EXON:           // NEG_TAG - NEG_AG
      pos1=begin; pos2=end; break;
    case NEG_SINGLE_EXON:          // NEG_TAG - NEG_ATG
      pos1=begin; pos2=end-3; break;
    }
}



BOOM::Vector<int> *GffPathFromParseGraph::getSignalCoordinates(
  BOOM::Vector<BOOM::GffTranscript*> &transcripts)
{
  BOOM::Vector<int> &path=*new BOOM::Vector<int>;
  int numTranscripts=transcripts.size();
  for(int i=0 ; i<numTranscripts ; ++i)
    {
      BOOM::GffTranscript &transcript=*transcripts[i];
      int numExons=transcript.getNumExons();
      for(int j=0 ; j<numExons ; ++j)
	{
	  BOOM::GffExon &exon=transcript.getIthExon(j);
	  Strand strand=exon.getStrand();
	  ContentType contentType=
	    ::exonTypeToContentType(exon.getExonType(),strand);
	  int begin=exon.getBegin(), end=exon.getEnd(), pos1, pos2;
	  exonCoordsToSigPos(contentType,begin,end,pos1,pos2);
	  path.push_back(pos1);
	  path.push_back(pos2);
	}
    }

  DirectComparator<int> cmp;
  BOOM::VectorSorter<int> sorter(path,cmp);
  sorter.sortAscendInPlace();
  return &path;
}


#endif
