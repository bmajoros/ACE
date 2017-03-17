/****************************************************************
 TranscriptSignals.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TranscriptSignals_H
#define INCL_TranscriptSignals_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/GffTranscript.H"
#include "BOOM/Essex.H"
#include "Signal.H"
#include "SignalSensors.H"
#include "StructureChange.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                    struct TranscriptSignal
 ****************************************************************/
struct TranscriptSignal {
  SignalType type; // GT, AG, TSS, TES
  float score, refScore, cutoff;
  int pos;
  bool broken;  // existed previously and was functional, but is now broken
  bool cryptic; // existed previously but wasn't annotated as functional
  bool denovo;  // didn't exist previously but does now
  bool weakened;// existed previously, is still functional but is weaker
  String seq;
  SignalType getType() { return type; }
  int getPos() { return pos; }
  void makeDead(); // a placeholder for exon-skipping/intron-retent event
  bool isDead();
  bool isBroken() { return broken; }
  bool isWeakened() { return weakened; }
  float getScore() { return score; }
  bool isCryptic() { return cryptic; }
  void setCryptic() { cryptic=true; }
  TranscriptSignal(SignalType,int pos,float score);
  TranscriptSignal() {}
  void printOn(ostream &);
};
ostream &operator<<(ostream &os,const TranscriptSignal &);



/****************************************************************
                     class TranscriptSignals
 ****************************************************************/
class TranscriptSignals {
public:
  TranscriptSignals();
  TranscriptSignals(const GffTranscript &);
  GffTranscript *toTranscript(const String &substrate,const SignalSensors &,
			      Essex::CompositeNode *&msg,bool reverrseStrand);
  int numSignals() const;
  TranscriptSignal &getIthSignal(int i);
  TranscriptSignal &operator[](int i) { return getIthSignal(i); }
  void deleteSignal(int index);
  TranscriptSignal &addSignal(SignalType,int pos,float score);
  void addSignal(const TranscriptSignal &);
  bool isCoding() const { return startCodon>-1; }
  void setStartCodon(int);
  int getStartCodon() const;
  bool anyBroken() const;
  void simulateBroken();
  bool anyWeakened() const;
  bool anyDead() const;
  bool anyCryptic() const;
  void printOn(ostream &) const;
  void setID(const String &i) { id=i; }
  void setSubstrate(const String &s) { substrate=s; }
  void setSource(const String &s) { source=s; }
  void setStrand(Strand s) { strand=s; }
  const String &getID() const { return id; }
  const String &getSubstrate() const { return substrate; }
  const String &getSource() const { return source; }
  Strand getStrand() const { return strand; }
  void setGeneID(const String &i) { geneID=i; }
  const String &getGeneID() const { return geneID; }
  StructureChange &getChange() { return structureChange; }
private:
  StructureChange structureChange;
  Vector<TranscriptSignal> signals;
  String id, geneID, substrate, source;
  Strand strand;
  int startCodon; // -1 if not defined (or noncoding gene)
};
ostream &operator<<(ostream &,const TranscriptSignals &);

#endif

