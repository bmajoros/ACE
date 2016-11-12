/****************************************************************
 EnumerateAltStructures.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "EnumerateAltStructures.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/Exceptions.H"
#include "BOOM/CombinationIterator.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                        AlternativeStructure
 ****************************************************************/

AlternativeStructure::AlternativeStructure(GffTranscript *t,ProteinFate fate)
  : transcript(t), proteinFate(fate), msg(NULL)
{
}



AlternativeStructure::~AlternativeStructure()
{
  delete transcript;
  delete msg;
}



void AlternativeStructure::reportCrypticSites(Essex::CompositeNode *parent,
					      bool reverseStrand,int L)
{
  for(Vector<TranscriptSignal>::const_iterator cur=crypticSignals.begin(),
	end=crypticSignals.end() ; cur!=end ; ++cur) {
    const TranscriptSignal &signal=*cur;
    Essex::CompositeNode *node=new Essex::CompositeNode("cryptic-site");
    String typeString=signal.getType()==GT ? "donor" : "acceptor";
    node->append(typeString);
    int pos=signal.getPos();
    if(reverseStrand) pos=L-pos-1;
    node->append(pos);
    node->append(signal.seq);
    node->append(signal.score);
    node->append("threshold:");
    node->append(signal.cutoff);
    parent->append(node);
  }  
}




/****************************************************************
                     EnumerateAltStructures
 ****************************************************************/

EnumerateAltStructures::EnumerateAltStructures(const TranscriptSignals &trans,
					       const String &genome,
					       int MAX_SPLICE_SHIFT,
					       int MIN_EXON_LEN,
					       int MIN_INTRON_LEN,
					       const int NMD_DISTANCE_PARM,
					       const SignalSensors &sensors,
					       bool allowExonSkipping,
					       bool allowIntronRetention,
					       bool allowCrypticSites)
  : original(trans), genome(genome), sensors(sensors),
    MAX_SPLICE_SHIFT(MAX_SPLICE_SHIFT), MIN_EXON_LEN(MIN_EXON_LEN),
    MIN_INTRON_LEN(MIN_INTRON_LEN), genomeSeq(genome,DnaAlphabet::global()),
    allowExonSkipping(allowExonSkipping), 
    allowIntronRetention(allowIntronRetention),
    allowCrypticSites(allowCrypticSites),
    nmd(NMD_DISTANCE_PARM)
{
  compute();
}



EnumerateAltStructures::~EnumerateAltStructures()
{
  for(Vector<AlternativeStructure*>::iterator cur=altStructures.begin(), 
	end=altStructures.end() ; cur!=end ; ++cur)
    delete *cur;
}



Vector<AlternativeStructure*> &EnumerateAltStructures::getAltStructures()
{
  return altStructures;
}



SignalSensor *EnumerateAltStructures::getSensor(SignalType type)
{
  switch(type) {
  case GT: return sensors.donorSensor;
  case AG: return sensors.acceptorSensor;
  default: INTERNAL_ERROR;
  }
}



void EnumerateAltStructures::findSites(SignalType type,int begin,int end,
				       const TranscriptSignal &notThisOne,
				       Vector<TranscriptSignal> &into)
{
  SignalSensor *sensor=getSensor(type);
  const int contextWindowLen=sensor->getContextWindowLength();
  const int consensusOffset=sensor->getConsensusOffset();
  for(int pos=begin ; pos<=end-contextWindowLen ; ++pos) {
    if(sensor->consensusOccursAt(genome,pos+consensusOffset)) {
      double score=sensor->getLogP(genomeSeq,genome,pos);
      if(score<sensor->getCutoff()) continue;
      TranscriptSignal signal(type,pos+consensusOffset,score);
      //if(signal.pos==notThisOne.pos) cout<<"skipping self..."<<endl;
      if(signal.pos==notThisOne.pos) continue;
      signal.setCryptic();
      /*signal.seq=genome.substring(pos,consensusOffset).tolower()+"_";
      signal.seq+=genome.substring(pos+consensusOffset,2)+"_";
      signal.seq+=genome.substring(pos+consensusOffset+2,
      contextWindowLen-consensusOffset-2).tolower();*/
      signal.cutoff=sensor->getCutoff();
      into.push_back(signal);
    }
  }
}



ContentType EnumerateAltStructures::classifyContent(TranscriptSignal prev,
						   TranscriptSignal next)
{
  switch(prev.getType()) {
  case TSS:
    switch(next.getType()) {
    case GT: return UTR5_INITIAL;
    case TES: return UTR5_SINGLE;
    default: INTERNAL_ERROR;
    } break;
  case GT: switch(next.getType()) {
    case AG: return UTR5_INTRON;
    default: INTERNAL_ERROR;
    }; break;
  case AG: switch(next.getType()) {
    case GT: return UTR5_INTERNAL;
    case TES: return UTR5_FINAL;
    default: INTERNAL_ERROR;
    } break;
  default: INTERNAL_ERROR;
  }
}



void EnumerateAltStructures::compute()
{
  // Generate lists of alternative splice sites for each broken one
  const int numSignals=original.numSignals();
  Array1D< Vector<TranscriptSignal> > alternatives(numSignals);
  for(int i=0 ; i<numSignals ; ++i) {
    TranscriptSignal signal=original[i];
    if(!signal.broken) { alternatives[i].push_back(signal); continue; }
    if(i==0 || i==numSignals-1) INTERNAL_ERROR; // TSS/TES can't be broken
    if(allowCrypticSites) {
      TranscriptSignal prev=original[i-1], next=original[i+1];
      const int prevSignalPos=prev.pos;
      const int nextSignalPos=next.pos;
      ContentType prevContent=classifyContent(prev,signal);
      ContentType nextContent=classifyContent(signal,next);
      int windowBegin=signal.pos-MAX_SPLICE_SHIFT;
      int windowEnd=signal.pos+MAX_SPLICE_SHIFT;
      const int leftLimit=prevSignalPos+
	(isIntron(prevContent) ? MIN_INTRON_LEN : MIN_EXON_LEN);
      const int rightLimit=nextSignalPos-
	(isIntron(nextContent) ? MIN_INTRON_LEN : MIN_EXON_LEN);
      if(windowBegin<leftLimit) windowBegin=leftLimit;
      if(windowEnd>rightLimit) windowEnd=rightLimit;
      if(windowBegin<windowEnd) {
	findSites(signal.getType(),windowBegin,windowEnd,signal,
		  alternatives[i]);
      }
    }
    if(allowExonSkipping || allowIntronRetention) {
      TranscriptSignal dead(signal.getType(),signal.getPos(),0.0);
      dead.makeDead();
      alternatives[i].push_back(dead);
    }
  }
  
  // Enumerate complete transcripts combinatorially
  CombinationIterator iter;
  for(int i=0 ; i<numSignals ; ++i) {
    const int n=alternatives[i].size();
    if(n<1) return;
    iter.addUnit(n);
  }
  for(iter.reset() ; iter.hasMore() ; iter.advance()) {
    Array1D<int> combination;
    iter.getCombination(combination);
    TranscriptSignals novel, temp;
    novel.setID(original.getID());
    novel.setGeneID(original.getGeneID());
    novel.setSubstrate(original.getSubstrate());
    novel.setSource(original.getSource());
    novel.setStrand(original.getStrand());
    novel.setStartCodon(original.getStartCodon());
    for(int i=0 ; i<numSignals ; ++i)
      novel.addSignal(alternatives[i][combination[i]]);
    bool deadSignals=novel.anyDead();
    if(deadSignals && allowIntronRetention) temp=novel;
    if(!deadSignals) addIfUnique(novel);
    else {
      if(allowExonSkipping && applyExonSkipping(novel))
	if(!novel.anyDead()) addIfUnique(novel);
      if(allowIntronRetention && deadSignals && applyIntronRetention(temp))
	if(!temp.anyDead()) addIfUnique(temp);
    }
  }
}


void EnumerateAltStructures::addIfUnique(TranscriptSignals signals)
{
  if(signals.anyCryptic()) signals.getChange().crypticSite=true;
  Essex::CompositeNode *msg=NULL;
  GffTranscript *transcript=signals.toTranscript(genome,sensors,msg);
  for(Vector<AlternativeStructure*>::iterator cur=altStructures.begin(), 
	end=altStructures.end() ; cur!=end ; ++cur) {
    const AlternativeStructure &other=**cur;
    if(other.transcript->identical(*transcript)) {
      delete transcript;
      return;
    }
  }
  int ejcDistance;
  ProteinFate fate=nmd.predict(*transcript,genome,ejcDistance);
  AlternativeStructure *structure=new AlternativeStructure(transcript,fate);
  structure->msg=msg;
  structure->ejcDistance=ejcDistance;
  structure->structureChange=signals.getChange();
  if(signals.anyCryptic())
    for(int i=0 ; i<signals.numSignals() ; ++i)
      if(signals[i].cryptic)  {
	structure->crypticSignals.push_back(signals[i]);
	//cout<<structure->crypticSignals.back().cutoff<<endl;
      }
  altStructures.push_back(structure);
}



bool EnumerateAltStructures::applyExonSkipping(TranscriptSignals &signals)
{
  bool changes=false;
  while(signals.anyDead()) {
    int numSignals=signals.numSignals();
    int i;
    for(i=2 ; i<numSignals-2 ; ++i) if(signals[i].isDead()) break;
    if(i>=numSignals-2) break;
    TranscriptSignal signal=signals[i];
    if(signal.getType()==GT)
      { signals.deleteSignal(i); signals.deleteSignal(i-1); changes=true; }
    else if(signal.getType()==AG)
      { signals.deleteSignal(i+1); signals.deleteSignal(i); changes=true; }
    else INTERNAL_ERROR;
  }
  if(changes) signals.getChange().exonSkipping=true;
  return changes;
}



bool EnumerateAltStructures::applyIntronRetention(TranscriptSignals &signals)
{
  bool changes=false;
  while(signals.anyDead()) {
    int numSignals=signals.numSignals();
    int i;
    for(i=0 ; i<numSignals ; ++i) if(signals[i].isDead()) break;
    if(i>=numSignals) INTERNAL_ERROR;
    TranscriptSignal signal=signals[i];
    if(signal.getType()==AG) 
      { signals.deleteSignal(i); signals.deleteSignal(i-1); changes=true; }
    else if(signal.getType()==GT)
      { signals.deleteSignal(i+1); signals.deleteSignal(i); changes=true; }
    else INTERNAL_ERROR;
  }
  if(changes) signals.getChange().intronRetention=true;
  return changes;
}









