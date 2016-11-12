/****************************************************************
 LightGraph.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_LightGraph_H
#define INCL_LightGraph_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/Regex.H"
#include "LightEdge.H"
#include "LightVertex.H"
using namespace std;
using namespace BOOM;

class LightGraph {
public:
  LightGraph(const String &filename);
  LightGraph(File &);
  virtual ~LightGraph();
  bool save(const String &filename);
  bool save(ostream &);
  int getNumVertices() const;
  int getNumEdges() const;
  int getLargestEdgeID() const;
  LightVertex *getVertex(int);
  LightEdge *getEdge(int);
  const String &getSubstrate() const;
  void printOn(ostream &);
  inline int getSubstrateLength() { return substrateLength; }
  void getATGs(Vector<LightVertex*> &);
  void getAnnotatedATGs(Vector<LightVertex*> &);
protected:
  Regex leftRegex, rightRegex, edgeIdRegex, annoRegex, IdRegex;
  Vector<LightVertex*> vertices;
  Vector<LightEdge*> edges;
  String substrate;
  int substrateLength;
  bool load(File &);
};

ostream &operator<<(ostream &,const LightGraph &);

#endif

