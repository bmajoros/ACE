/****************************************************************
 TranscriptPaths.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "TranscriptPaths.H"
#include "NBest.H"
#include "BOOM/Stack.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/VectorSorter.H"
using namespace std;
using namespace BOOM;

TranscriptPaths::TranscriptPaths(LightGraph &G,int maxPaths,int seqLen,
				 const Model &model)
  : G(G), seqLen(seqLen), model(model)
{
  buildPaths(maxPaths);
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



void TranscriptPaths::buildPaths(int N)
{
  // Extract N best paths using dynamic programming
  //cout<<"running N-best"<<endl;
  NBest nbest(G,N);
  //cout<<"N-best traceback"<<endl;
  nbest.getPaths(paths);

  // Score the paths
  //cout<<"scoring paths "<<paths.size()<<endl;
  for(Vector<TranscriptPath*>::iterator cur=paths.begin(), end=paths.end() ;
      cur!=end ; ++cur) 
    (*cur)->computeScore(model);
}



void TranscriptPaths::computeLRs(double denom)
{
  const double L=double(seqLen);
  int numPaths=paths.size();
  for(int i=0 ; i<numPaths ; ++i) {
    TranscriptPath *path=paths[i];
    path->computeScore(model);
    const double llr=path->getScore()/L-denom;
    path->setScore(exp(llr));
  }
}



// ### Not currently used:
double sigmoid(double x)
{
  return 1.0/(1.0+exp(-x));
}



void TranscriptPaths::computePosteriors()
{
  // First, collect log probabilities, log(P(path,seq))
  Vector<double> logProbs;
  int numPaths=paths.size();
  for(int i=0 ; i<numPaths ; ++i) {
    TranscriptPath *path=paths[i];
    path->computeScore(model);
    logProbs.push_back(path->getScore());}
  
  // Marginalize out the paths to get log(P(seq))
  double logSum=logProbs.empty() ? 0.0 : sumLogProbs(logProbs);
  if(!isFinite(logSum)) {
    cerr<<"logSum="<<logSum<<endl;
    for(int i=0 ; i<numPaths ; ++i) cerr<<" score="<<paths[i]->getScore();
    cerr<<endl;
    INTERNAL_ERROR;}

  // Convert P(path,seq) to P(path|seq) by dividing by P(seq)
  double sum=0.0;
  for(int i=0 ; i<numPaths ; ++i) {
    TranscriptPath *path=paths[i];
    double newScore=exp(path->getScore()-logSum);
    //double newScore=path->getScore()-logSum;
    //double newScore=sigmoid(path->getScore());
    if(!isFinite(newScore)) {
      cerr<<"newScore="<<newScore<<endl;
      INTERNAL_ERROR;}
    path->setScore(newScore);
    sum+=newScore;
  }

  // Normalize again in non-log space, to fix rounding issues
  if(false) {
    for(int i=0 ; i<numPaths ; ++i) {
      TranscriptPath *path=paths[i];
      path->setScore(path->getScore()/sum);
      if(!isFinite(path->getScore())) {
	cerr<<"path score="<<path->getScore()<<endl;
	INTERNAL_ERROR;}}}
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

