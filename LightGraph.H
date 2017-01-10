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
  LightGraph(const String &substrate,int substrateLength);
  virtual ~LightGraph();
  void addVertex(LightVertex*);
  void addEdge(LightEdge*);
  void sort();
  void sortVertices();
  void sortEdges();
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
  void dropVertex(int index); // doesn't delete it, sets pointer to NULL
  void dropEdge(int index); // doesn't delete it, sets pointer to NULL
  void deleteNullVertices();
  void deleteNullEdges();
  void deleteDuplicates();
  void deleteVertex(int index); // dangerous! invalidates indices!
  void deleteEdge(int index); // dangerous! invalidates indices!
  bool vertexExists(const String &substrate,Strand,int begin,int end,
		    SignalType) const;
  bool edgeExists(const String &substrate,Strand,int begin,int end,
		  ContentType) const;
protected:
  Regex leftRegex, rightRegex, edgeIdRegex, annoRegex, IdRegex;
  Vector<LightVertex*> vertices;
  Vector<LightEdge*> edges;
  String substrate;
  int substrateLength;
  bool load(File &);
  void deleteIncidentEdges(LightVertex *);
  void deleteEdges(Vector<LightEdge*> &);
};

ostream &operator<<(ostream &,const LightGraph &);

#endif

