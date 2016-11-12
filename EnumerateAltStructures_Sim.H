/****************************************************************
 EnumerateAltStructures.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_EnumerateAltStructures_H
#define INCL_EnumerateAltStructures_H
#include <iostream>
#include "BOOM/GffTranscript.H"
#include "BOOM/Vector.H"
#include "BOOM/Sequence.H"
#include "BOOM/Essex.H"
#include "NMD.H"
#include "TranscriptSignals.H"
#include "Isochore.H"
#include "SignalSensors.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                     struct AlternativeStructure
 ****************************************************************/
struct AlternativeStructure {
  GffTranscript *transcript;
  ProteinFate proteinFate;
  int ejcDistance;
  StructureChange structureChange;
  Essex::CompositeNode *msg;
  Vector<TranscriptSignal> crypticSignals;
  AlternativeStructure(GffTranscript *,ProteinFate);
  ~AlternativeStructure();
  void reportCrypticSites(Essex::CompositeNode *,bool reverseStrand,
			  int substrateLen);
};



/****************************************************************
                     class EnumerateAltStructures
 ****************************************************************/
class EnumerateAltStructures {
public:
  EnumerateAltStructures(const TranscriptSignals &original,
			 const String &altGenome,int MAX_SPLICE_SHIFT,
			 int MIN_EXON_LEN,int MIN_INTRON_LEN,
			 const int NMD_DISTANCE_PARM,
			 const SignalSensors &sensors,
			 bool allowExonSkipping,bool allowIntronRetention,
			 bool allowCrypticSites);
  virtual ~EnumerateAltStructures();
  Vector<AlternativeStructure*> &getAltStructures();
private:
  Vector<AlternativeStructure*> altStructures;
  int MAX_SPLICE_SHIFT, MIN_EXON_LEN, MIN_INTRON_LEN;
  const TranscriptSignals &original;
  const String &genome;
  const Sequence genomeSeq;
  const SignalSensors &sensors;
  NMD nmd;
  bool allowExonSkipping, allowIntronRetention, allowCrypticSites;
  void compute();
  void findSites(SignalType,int begin,int end,const TranscriptSignal 
		 &notThisOne,Vector<TranscriptSignal> &into);
  SignalSensor *getSensor(SignalType);
  ContentType classifyContent(TranscriptSignal prev,TranscriptSignal next);
  bool applyExonSkipping(TranscriptSignals &);
  bool applyIntronRetention(TranscriptSignals &);
  void addIfUnique(TranscriptSignals);
};

#endif

