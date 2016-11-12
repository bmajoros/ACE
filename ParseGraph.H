/****************************************************************
 ParseGraph.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ParseGraph_H
#define INCL_ParseGraph_H
#include <iostream>
#include "BOOM/Stack.H"
#include "BOOM/Array1D.H"
#include "Signal.H"
#include "genezilla.H"

class GeneZilla;

#ifdef EXPLICIT_GRAPHS
class ParseGraph
{
  BOOM::Vector<SignalPtr> *vertices;
  GeneZilla &genezilla;
public:
  ParseGraph(GeneZilla &);
  virtual ~ParseGraph();
  void addVertex(SignalPtr);
  int numVertices();
  SignalPtr getIthVertex(int);
  BOOM::Stack<SignalPtr> *findOptimalPath(double &parseScore);
  void deleteVertices(BOOM::Array1D<bool> &deleteIfTrue);
  void setVertexIndices();
  int findSignal(SignalType,int consensusPos); // -1 if not found
  int findFirstSignalAfter(int pos); // -1 if no signals follow pos
  int countLeftTermini();
  int countRightTermini();
  void toGFF(ostream &) const;
  int countEdges() const;
  //void normalize();
};

ostream &operator<<(ostream &,const ParseGraph &);

#endif

#endif
