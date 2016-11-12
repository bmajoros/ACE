/****************************************************************
 TranscriptSignals.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/Constants.H"
#include "TranscriptSignals.H"
#include "StartCodonFinder.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                    struct TranscriptSignal
 ****************************************************************/
TranscriptSignal::TranscriptSignal(SignalType t,int pos,float score)
  : type(t), score(score), pos(pos), broken(false), cryptic(false)
{
  // ctor
}



void TranscriptSignal::makeDead()
{
  score=NEGATIVE_INFINITY;
}



bool TranscriptSignal::isDead()
{
  return score==NEGATIVE_INFINITY;
}



void TranscriptSignal::printOn(ostream &os)
{
  os<<type<<":"<<pos<<":"<<score;
}



ostream &operator<<(ostream &os,const TranscriptSignal &sig)
{
  sig.printOn(os);
  return os;
}



/****************************************************************
                     class TranscriptSignals
 ****************************************************************/
TranscriptSignals::TranscriptSignals()
  : strand(FORWARD_STRAND), startCodon(-1)
{
  // default ctor
}



TranscriptSignals::TranscriptSignals(const GffTranscript &transcript)
{
  id=transcript.getTranscriptId();
  source=transcript.getSource();
  substrate=transcript.getSubstrate();
  geneID=transcript.getGeneId();
  if(transcript.getStrand()!=FORWARD_STRAND) 
    throw "please reverse-complement the sequence first";
  strand=FORWARD_STRAND;
  Vector<GffExon*> exons;
  transcript.getRawExons(exons);
  const int numExons=exons.size();
  for(int i=0 ; i<numExons ; ++i) {
    const GffExon *exon=exons[i];
    int exonBegin=exon->getBegin(), exonEnd=exon->getEnd();
    if(i==0) { // initial exon
      addSignal(TSS,exonBegin,0.0);
      addSignal(GT,exonEnd,0.0);
    }
    else if(i==numExons-1) { // final exon
      addSignal(AG,exonBegin-2,0.0);
      addSignal(TES,exonEnd,0.0);
    }
    else { // internal exon
      addSignal(AG,exonBegin-2,0.0);
      addSignal(GT,exonEnd,0.0);
    }
  }
  for(Vector<GffExon*>::iterator cur=exons.begin(), end=exons.end() ;
      cur!=end ; ++cur)  delete *cur;
}



void TranscriptSignals::simulateBroken()
{
  const int numSignals=signals.size();
  if(numSignals<=2) return;
  const int index=1+RandomNumber(numSignals-2);
  TranscriptSignal &signal=signals[index];
  signal.broken=true;
}



GffTranscript *TranscriptSignals::toTranscript(const String &genome,
					       const SignalSensors &sensors,
					       Essex::CompositeNode *&msg,
					       bool reverseStrand)
{
  msg=NULL;
  GffTranscript *transcript=new GffTranscript(id,substrate,strand,source);
  transcript->setGeneId(geneID);
  const int numSignals=signals.size();
  if(numSignals%2) INTERNAL_ERROR; // numSignals must be even
  if(isCoding()) { // protein-coding gene
    for(int i=0 ; i<numSignals ; i+=2) {
      TranscriptSignal thisSig=signals[i], nextSig=signals[i+1];
      int begin=thisSig.pos, end=nextSig.pos;
      if(thisSig.type==AG) begin+=2;
      GffExon *exon=
	new GffExon(ET_EXON,begin,end,*transcript,false,0.0,false,0);
      transcript->addExon(exon);
    }
    int newStart=StartCodonFinder::findStartCodon(*transcript,
						  transcript->peekExons(),
						  genome,
						  startCodon,
						  sensors);
    if(newStart!=startCodon)
      if(newStart>0) {
	msg=new Essex::CompositeNode("start-codon-change");
	if(reverseStrand) {
	  msg->append("from",int(genome.length()-startCodon));
	  msg->append("to",int(genome.length()-newStart));}
	else {
	  msg->append("from",startCodon);
	  msg->append("to",newStart);}
      }
      else msg=new Essex::CompositeNode("start-codon-lost");
    if(newStart>=0)
      transcript->splitUTRandCDS(genome,newStart,sensors.stopCodons);
  }
  else { // noncoding gene
    for(int i=0 ; i<numSignals ; i+=2) {
      TranscriptSignal thisSig=signals[i], nextSig=signals[i+1];
      int begin=thisSig.pos, end=nextSig.pos;
      if(thisSig.type==AG) begin+=2;
      GffExon *exon=
	new GffExon(ET_UTR5,begin,end,*transcript,false,0.0,false,0);
      transcript->addUTR(exon);
    }
  }
  return transcript;
}



int TranscriptSignals::numSignals() const
{
  return signals.size();
}



TranscriptSignal &TranscriptSignals::getIthSignal(int i)
{
  return signals[i];
}



TranscriptSignal &TranscriptSignals::addSignal(SignalType t,int pos,
					       float score)
{
  signals.push_back(TranscriptSignal(t,pos,score));
  return signals.back();
}



void TranscriptSignals::addSignal(const TranscriptSignal &signal)
{
  signals.push_back(signal);
}



void TranscriptSignals::setStartCodon(int s)
{
  startCodon=s;
}



int TranscriptSignals::getStartCodon() const
{
  return startCodon;
}



bool TranscriptSignals::anyBroken() const
{
  for(Vector<TranscriptSignal>::const_iterator cur=signals.begin(), end=
	signals.end() ; cur!=end ; ++cur)
    if((*cur).isBroken()) return true;
  return false;
}



bool TranscriptSignals::anyWeakened() const
{
  for(Vector<TranscriptSignal>::const_iterator cur=signals.begin(), end=
	signals.end() ; cur!=end ; ++cur)
    if((*cur).isWeakened()) return true;
  return false;
}



bool TranscriptSignals::anyCryptic() const
{
  for(Vector<TranscriptSignal>::const_iterator cur=signals.begin(), end=
	signals.end() ; cur!=end ; ++cur)
    if((*cur).isCryptic()) return true;
  return false;
}



bool TranscriptSignals::anyDead() const
{
  for(Vector<TranscriptSignal>::const_iterator cur=signals.begin(), end=
	signals.end() ; cur!=end ; ++cur)
    if((*cur).isDead()) return true;
  return false;
}



void TranscriptSignals::deleteSignal(int index)
{
  signals.cut(index);
}



void TranscriptSignals::printOn(ostream &os) const
{
  os<<"[start:"<<startCodon;
  for(Vector<TranscriptSignal>::const_iterator cur=signals.begin(), end=
	signals.end() ; cur!=end ; ++cur)
    os<<"|"<<*cur;
  os<<"]";
}



ostream &operator<<(ostream &os,const TranscriptSignals &signals)
{
  signals.printOn(os);
  return os;
}




