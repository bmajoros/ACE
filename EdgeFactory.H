/****************************************************************
 EdgeFactory.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_EdgeFactory_H
#define INCL_EdgeFactory_H
#include "Edge.H"
#include "EvidenceFilter.H"


class EdgeFactory {
public:
  virtual ~EdgeFactory() {}
  virtual PhasedEdge *newPhasedEdge(double scorePhase0,double scorePhase1,
				    double scorePhase2,SignalPtr left,
				    SignalPtr right);
  virtual NonPhasedEdge *newNonPhasedEdge(double score,SignalPtr left,
					  SignalPtr right);
};


class FilteredEdgeFactory : public EdgeFactory {
public:
  FilteredEdgeFactory(EvidenceFilter *);
  virtual ~FilteredEdgeFactory();
  virtual PhasedEdge *newPhasedEdge(double scorePhase0,double scorePhase1,
				    double scorePhase2,SignalPtr left,
				    SignalPtr right);
  void setEvidence(EvidenceFilter *);
protected:
  EvidenceFilter *filter;
};


#endif
