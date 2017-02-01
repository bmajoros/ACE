/****************************************************************
 NBest.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/Constants.H"
#include "BOOM/Stack.H"
#include "NBest.H"
using namespace std;
using namespace BOOM;



NBest::NBest(LightGraph &G,int N)
  : G(G), N(N)
{
  buildTrellis();
}



void NBest::buildTrellis()
{
  int queueCapacity=N<100 ? N : 100;
  const int numVertices=G.getNumVertices();
  links.resize(numVertices);
  if(numVertices==0) return;
  for(int i=0 ; i<numVertices ; ++i) {
    FixedSizePriorityQueue<TrellisLink> Q(queueCapacity,cmp);
    LightVertex *signal=G.getVertex(i);
    if(!signal) continue;
    Vector<LightEdge*> &edges=signal->getEdgesIn();
    int numEdges=edges.size();
    if(numEdges==0) Q.insert(TrellisLink(NULL,NULL)); // left terminus
    for(Vector<LightEdge*>::iterator cur=edges.begin(), end=edges.end() ; 
	cur!=end ; ++cur) {
      LightEdge *currentEdge=*cur;
      LightVertex *pred=currentEdge->getLeft();
      Array1D<TrellisLink> &predLinks=links[pred->getID()];
      const int numPredLinks=predLinks.size();
      for(int j=0 ; j<numPredLinks ; ++j) {
	TrellisLink &predLink=predLinks[j];
	const double score=currentEdge->getScore()+predLink.getScore();
	if(isFinite(score))
	  Q.insert(TrellisLink(&predLink,currentEdge,score)); } }

    // Copy selected links into signal's link set
    Array1D<TrellisLink> &linkSet=links[i];
    const int size=Q.getNumElements();
    linkSet.resize(size);
    int j=0;
    for(FixedSizePriorityQueue<TrellisLink>::iterator cur=Q.begin(), 
	  end=Q.end() ; cur!=end ; ++cur) linkSet[j++]=*cur;
  }
  
  // Find all of the right-terminal links
  for(int rCur=numVertices-1; rCur>=0; --rCur) {
    LightVertex *rt=G.getVertex(rCur);
    if(!rt) continue;
    if(!rt->getEdgesOut().isEmpty()) break;
    Array1D<TrellisLink> &rtLinks=links[rCur];
    const int numRtLinks=rtLinks.size();
    for(int j=0 ; j<numRtLinks ; ++j)
      termini.push_back(rtLinks[j]);
  }
  
  // Return links representing the N best parses
  VectorSorter<TrellisLink> sorter(termini,cmp);
  sorter.sortDescendInPlace();
  if(N<termini.size()) termini.resize(N);
}



void NBest::getPaths(Vector<TranscriptPath*> &into)
{
  const int numTermini=termini.size();
  for(int i=0 ; i<numTermini ; ++i) {
    TranscriptPath *path=new TranscriptPath();
    traceback(&termini[i],*path);
    into.push_back(path);
  }
}



void NBest::traceback(TrellisLink *endLink,TranscriptPath &path) 
{
  Stack<TrellisLink*> pStack;
  TrellisLink *currentLink=endLink;
  while(currentLink) {
    pStack.push(currentLink);
    currentLink=currentLink->getPred();
  }
  while(!pStack.isEmpty()) {
    TrellisLink *pl=pStack.pop();
    path.addEdge(dynamic_cast<ACEplus_Edge*>(pl->getEdge()));
  }
}
