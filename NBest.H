/****************************************************************
 NBest.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_NBest_H
#define INCL_NBest_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/FixedSizePriorityQueue.H"
#include "BOOM/Array1D.H"
#include "LightGraph.H"
#include "TranscriptPath.H"
#include "TrellisLink.H"
using namespace std;
using namespace BOOM;

/****************************************************************
 class NBest
 ****************************************************************/
class NBest {
public:
  NBest(LightGraph &,int N);
  void getPaths(Vector<TranscriptPath*> &into);
private:
  LightGraph &G;
  const int N;
  Array1D< Array1D<TrellisLink> > links;
  BOOM::Vector<TrellisLink> termini;
  TrellisLinkComparator cmp;
  void buildTrellis();
  void traceback(TrellisLink *endLink,TranscriptPath &into);
};

#endif

