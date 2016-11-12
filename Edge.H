/****************************************************************
 Edge.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_Edge_H
#define INCL_Edge_H
#include "BOOM/Interval.H"
#include "Signal.H"
#include "genezilla.H"
#include "ContentType.H"
#include "BOOM/Strand.H"

#define SignalPtr BOOM::SmartPointer<Signal>

//=================================================================
//                           class Edge
//=================================================================

/*
  Abstract base class.
*/

class Edge
{
protected:
  SignalPtr left, right;
  Edge(SignalPtr left,SignalPtr right);
public:
  virtual ~Edge();
  SignalPtr getLeft();
  SignalPtr getRight();
  virtual void dropScores()=0;
  virtual bool isPhased() const=0;
  int getEdgeBegin(); // after last base of left signal's context window
  int getEdgeEnd();   // first base of right signal's context window
  int getFeatureBegin(); // specific to exon/intron/UTR/intergenic
  int getFeatureEnd();   // specific to exon/intron/UTR/intergenic
  int getFeatureLength() {return getFeatureEnd()-getFeatureBegin();}
  Interval getFeatureInterval() { return Interval(getFeatureBegin(),getFeatureEnd()); }
  ContentType getContentType() const;
  Strand getStrand();
  bool isCoding();
  bool isIntron();
  bool isIntergenic();
  bool isUTR();
  int propagateForward(int phase);
  int propagateBackward(int phase);
  virtual void printOn(ostream &) const;
};
ostream &operator<<(ostream &,const Edge &);


//=================================================================
//                        class PhasedEdge
//=================================================================

/*
  Used for exons and introns.  Note that the phase of an edge is the
  phase of the preceding signal (i.e., its left vertex).
 */

class PhasedEdge : public Edge
{
  double scores[3]; // these include duration and transition scores
public:
  PhasedEdge(double scorePhase0,double scorePhase1,double scorePhase2,
	     SignalPtr left,SignalPtr right);
  double getEdgeScore(int phase); // phase is of left signal
  void setEdgeScore(int phase,double);
  double getInductiveScore(int phase); // phase is of left signal
  virtual bool isPhased() const;
  virtual void printOn(ostream &) const;
  virtual void dropScores() { for(int i=0 ; i<3 ; ++i) scores[i]=0; }
};



//=================================================================
//                       class NonPhasedEdge
//=================================================================

/*
  Used for intergenic & UTR
*/

class NonPhasedEdge : public Edge
{
  double score; // includes duration and transition scores
public:
  NonPhasedEdge(double score,SignalPtr left,SignalPtr right);
  double getEdgeScore();
  void setEdgeScore(double);
  double getInductiveScore();
  virtual bool isPhased() const;
  virtual void printOn(ostream &) const;
  virtual void dropScores() { score=0; }
};


#endif
