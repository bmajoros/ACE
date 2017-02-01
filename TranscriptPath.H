/****************************************************************
 TranscriptPath.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TranscriptPath_H
#define INCL_TranscriptPath_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/GffTranscript.H"
#include "BOOM/Comparator.H"
#include "ACEplus_Edge.H"
#include "ACEplus_Vertex.H"
using namespace std;
using namespace BOOM;


class TranscriptPath {
public:
  TranscriptPath();
  void addEdge(ACEplus_Edge*);
  int numEdges() const;
  ACEplus_Edge *operator[](int);
  ACEplus_Edge *lastEdge();
  ACEplus_Vertex *lastVertex();
  double getScore() const { return score; }
  void computeScore();
  void setScore(double s) { score=s; }
  TranscriptPath *clone();
  const StructureChange &getChange() const;
  GffTranscript *toTranscript(const String &transcriptID,
			      const String &geneID,
			      const String &substrate,
			      Strand,
			      const String &source) const;
  void getVertices(Vector<ACEplus_Vertex*> &into);
  void dumpScores();
  bool isFullyAnnotated() const;
protected:
  Vector<ACEplus_Edge*> edges;
  double score;
  StructureChange change;
};

class TranscriptPathComparator : public Comparator<TranscriptPath*> {
public:
  virtual bool equal(TranscriptPath *&a,TranscriptPath *&b) 
    { return a->getScore()==b->getScore(); }
  virtual bool greater(TranscriptPath *&a,TranscriptPath *&b) 
    { return a->getScore()>b->getScore(); }
  virtual bool less(TranscriptPath *&a,TranscriptPath *&b) 
    { return a->getScore()<b->getScore(); }
};

#endif

