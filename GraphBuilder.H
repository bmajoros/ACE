/****************************************************************
 GraphBuilder.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_GraphBuilder_H
#define INCL_GraphBuilder_H
#include <iostream>
#include "BOOM/GffTranscript.H"
#include "BOOM/String.H"
#include "BOOM/Sequence.H"
#include "BOOM/Array1D.H"
#include "BOOM/CigarAlignment.H"
#include "BOOM/Direction.H"
#include "BOOM/Set.H"
#include "TranscriptSignals.H"
#include "LightGraph.H"
#include "Model.H"
#include "ACEplus_Vertex.H"
#include "ACEplus_Edge.H"
using namespace std;
using namespace BOOM;

/* This class assumes genes have been mapped to the forward strand */

class GraphBuilder {
public:
  GraphBuilder(const GffTranscript &projected,
	       const TranscriptSignals &,
	       Model &,
	       const Sequence &refSeq,const String &refSeqStr,
	       const Sequence &altSeq,const String &altSeqStr,
	       CigarAlignment &altToRef,bool strict=false);
  virtual ~GraphBuilder() {}
  LightGraph *getGraph();
  bool mapped() const;
protected:
  struct ExonEdge {
    LightEdge *edge;
    Interval interval;
  };
  const GffTranscript &projected;
  const TranscriptSignals &signals;
  const Sequence &refSeq, &altSeq;
  const String &refSeqStr, &altSeqStr;
  const CigarAlignment &altToRef;
  bool changes;
  Model &model;
  LightGraph *G;
  Vector<Interval> variants;
  bool buildGraph(bool strict);
  bool buildTranscriptGraph();
  ACEplus_Vertex *newVertex(const String &substrate,SignalType,int begin,
			    int end,double score,Strand,int ID,
			    bool denovo=false);
  ACEplus_Edge *newEdge(const String &substrate,ContentType,LightVertex *,
		     LightVertex *,int begin,int end,Strand,int ID);
  double scoreSignal(SignalType,int pos,const Sequence &,const String &);
  double scoreEdge(ACEplus_Edge *);
  void getSignalWindow(SignalType,int consensusPos,int &windowBegin,
		       int &windowEnd,int &consensusLength);
  ContentType getContentType(SignalType,SignalType);
  void handleBrokenSites();
  void handleBrokenSite(LightVertex *);
  void handleBrokenSite_cryptic(LightVertex *);
  void handleBrokenSite_skipping(LightVertex *);
  void handleBrokenSite_retention(LightVertex *);
  void handleSkippingLeft(LightVertex *);
  void handleSkippingRight(LightVertex *);
  void handleRetentionLeft(LightVertex *);
  void handleRetentionRight(LightVertex *);
  void leftSweep(LightGraph &G,Array1D<int> &vertexCounts,
		 Array1D<int> &edgeCounts);
  void rightSweep(LightGraph &G,Array1D<int> &vertexCounts,
		  Array1D<int> &edgeCounts);
  void pruneUnreachable(LightGraph &G);
  void deleteUnreachable(LightGraph &G,
			 Array1D<int> &vertexLeftCounts,
			 Array1D<int> &vertexRightCounts,
			 Array1D<int> &edgeLeftCounts,
			 Array1D<int> &edgeRightCounts);
  void getContextWindow(LightVertex *v,
			int &begin,
			int &end);
  void handleDeNovoSites();
  void findVariants(Vector<Interval> &variants);
  void scanDeNovo(SignalSensor &,Vector<Interval> &variants,
		  Vector<ACEplus_Vertex*> &newVertices);
  void getVariantWindows(Vector<Interval> &windows,SignalSensor &,
			 Vector<Interval> &variants);
  void coalesceWindows(Vector<Interval> &windows);
  void deNovoSignalSensing(SignalSensor &,
			   Vector<Interval> &windows,
			   Vector<ACEplus_Vertex*> &newVertices);
  void linkDeNovoVertices(Vector<ACEplus_Vertex*> &newVertices);
  void linkDeNovoLeft(ACEplus_Vertex *,int id);
  void linkDeNovoRight(ACEplus_Vertex *,int id);
  void getLinkTypes(SignalType,Direction,Set<SignalType> &linkTypes);
  void addCrypticExons(Vector<ACEplus_Vertex*> &newVertices);
  bool allVerticesAreAnnotated();
  void scanCrypticExonLeft(int of);
  void scanCrypticExonRight(int of);
  LightVertex *findAnnotatedVertexLeftOf(int of);
  LightVertex *findAnnotatedVertexRightOf(int of);
  void findMateSignals(SignalSensor &,const Interval &scanWindow,
		       Vector<ACEplus_Vertex*> &newVertices);
  ACEplus_Edge *linkVertices(LightVertex *left,LightVertex *right);
  void handleExonWeakening(const Vector<Interval> &,Vector<ExonEdge> &exons);
  void handleExonExtension(const Vector<Interval> &,Vector<ExonEdge> &exons);
  void handleExonShortening(const Vector<Interval> &variants,
			    Vector<ExonEdge> &exons);
  void handleDeepIntronic(const Vector<Interval> &,Vector<ExonEdge> &exons);
  void handleRegulatoryChanges();
  double exonDefChange(Interval); // alt vs. ref
  double exonIntronRatio(Interval interval); // exon vs. intron
  double exonIntronRatio(LightVertex *from,LightVertex *to);
  void getAnnotatedExons(Vector<ExonEdge> &into);
  void getVariantsInExons(const Vector<Interval> &variants,
			  const Vector<ExonEdge> &exons,
			  const Vector< pair<Interval,ExonEdge> > &into);
  void addExonSkippingEdge(ExonEdge exon);
  void getAllVertices(int fromID,int toID,SignalType,Vector<int> &into);
  void lengthenExon(const Interval &variant,const ExonEdge &);
  void lengthenExonLeft(const Interval &variant,const ExonEdge &);
  void lengthenExonRight(const Interval &variant,const ExonEdge &);
  void scan(const Interval &region,SignalType,Vector<ACEplus_Vertex*> &into);
  bool isInExon(Interval variant,Vector<ExonEdge> &exons);
  int distanceOfExonicVariantToEnd(const Interval &variant,const ExonEdge &);
  void shortenExon(const Interval &variant,const ExonEdge &exon);
  void shortenExonFromTheRight(const Interval &variant,const ExonEdge &exon);
  void shortenExonFromTheLeft(const Interval &variant,const ExonEdge &exon);
  void proposeCrypticExons(const Interval &variant,
			   const ExonEdge &prevExon,
			   const ExonEdge &nextExon);
  void inspectEdgeScore(ACEplus_Edge *,const String &label);
  LightEdge *getAnnotatedEdge(Vector<LightEdge*> &);
  void getAnnotatedEdges(Vector<LightEdge*> &into);
  bool coversAnnotatedIntron(const LightEdge &edge,
			     Vector<LightEdge*> &annotatedEdges);
  void markIntronRetentions();
};

#endif

