/****************************************************************
 TranscriptPaths.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "TranscriptPaths.H"
#include "BOOM/Stack.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/VectorSorter.H"
using namespace std;
using namespace BOOM;

TranscriptPaths::TranscriptPaths(LightGraph &G)
  : G(G)
{
  buildPaths();
}


TranscriptPaths::~TranscriptPaths()
{
  for(Vector<TranscriptPath*>::iterator cur=paths.begin(), end=paths.end() ;
      cur!=end ; ++cur) delete *cur;
}



int TranscriptPaths::numPaths() const
{
  return paths.size();
}



TranscriptPath *TranscriptPaths::operator[](int i)
{
  return paths[i];
}



void TranscriptPaths::buildPaths()
{
  const int numVertices=G.getNumVertices();
  if(numVertices<2) return;

  // Initialize recursion stack
  Stack<TranscriptPath*> S;
  LightVertex *leftTerminus=G.getVertex(0);
  LightVertex *rightTerminus=G.getVertex(numVertices-1);
  Vector<LightEdge*> &out=leftTerminus->getEdgesOut();
  for(Vector<LightEdge*>::iterator cur=out.begin(), end=out.end() ;
      cur!=end ; ++cur) {
    LightEdge *edge=*cur;
    TranscriptPath *path=new TranscriptPath();
    path->addEdge(dynamic_cast<ACEplus_Edge*>(edge));
    S.push(path);
  }

  // Perform depth-first search using stack
  while(!S.isEmpty()) {
    TranscriptPath *path=S.pop();
    LightVertex *v=path->lastVertex();
    if(v==rightTerminus) { paths.push_back(path); continue; }
    Vector<LightEdge*> &out=v->getEdgesOut();
    for(Vector<LightEdge*>::iterator cur=out.begin(), end=out.end() ;
	cur!=end ; ++cur) {
      LightEdge *edge=*cur;
      TranscriptPath *newPath=path->clone();
      newPath->addEdge(dynamic_cast<ACEplus_Edge*>(edge));
      S.push(newPath);
    }
    delete path;
  }

  // Score the paths
  for(Vector<TranscriptPath*>::iterator cur=paths.begin(), end=paths.end() ;
      cur!=end ; ++cur) 
    (*cur)->computeScore();
}



void TranscriptPaths::computeLRs(double denom)
{
  int numPaths=paths.size();
  for(int i=0 ; i<numPaths ; ++i) {
    TranscriptPath *path=paths[i];
    path->computeScore();
    //cout<<"path->getScore()="<<path->getScore()<<endl;
    const double llr=path->getScore()-denom;
    const double lr=exp(llr);
    if(lr>100000)
      cout<<"setting score "<<lr<<" = "<<path->getScore()<<" - "<<denom<<endl;
    path->setScore(lr);
  }
}



void TranscriptPaths::computePosteriors()
{
  // First, collect log probabilities, log(P(path,seq))
  Vector<double> logProbs;
  int numPaths=paths.size();
  for(int i=0 ; i<numPaths ; ++i) {
    TranscriptPath *path=paths[i];
    path->computeScore();
    logProbs.push_back(path->getScore());
  }
  
  // Marginalize out the paths to get log(P(seq))
  double logSum=sumLogProbs(logProbs);
  if(!isFinite(logSum)) {
    cerr<<"logSum="<<logSum<<endl;
    for(int i=0 ; i<numPaths ; ++i) cerr<<" score="<<paths[i]->getScore();
    cerr<<endl;
    INTERNAL_ERROR;
  }

  // Convert P(path,seq) to P(path|seq) by dividing by P(seq)
  double sum=0.0;
  for(int i=0 ; i<numPaths ; ++i) {
    TranscriptPath *path=paths[i];
    double newScore=exp(path->getScore()-logSum);
    if(!isFinite(newScore)) {
      cerr<<"newScore="<<newScore<<endl;
      INTERNAL_ERROR;
    }
    path->setScore(newScore);
    sum+=newScore;
  }

  // Normalize again in non-log space, to fix rounding issues
  for(int i=0 ; i<numPaths ; ++i) {
    TranscriptPath *path=paths[i];
    path->setScore(path->getScore()/sum);
    if(!isFinite(path->getScore())) {
      cerr<<"path score="<<path->getScore()<<endl;
      INTERNAL_ERROR;
    }
  }
}



void TranscriptPaths::filter(double minScore)
{
  int n=paths.size();
  for(int i=0 ; i<n ; ++i)
    if(paths[i]->getScore()<minScore) {
      paths.cut(i);
      --i;
      --n;
    }
}



void TranscriptPaths::sort()
{
  TranscriptPathComparator cmp;
  VectorSorter<TranscriptPath*> sorter(paths,cmp);
  sorter.sortDescendInPlace();
}

