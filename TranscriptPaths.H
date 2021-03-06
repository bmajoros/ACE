/****************************************************************
 TranscriptPaths.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TranscriptPaths_H
#define INCL_TranscriptPaths_H
#include <iostream>
#include "BOOM/Vector.H"
#include "LightGraph.H"
#include "TranscriptPath.H"
#include "Model.H"
using namespace std;
using namespace BOOM;


class TranscriptPaths {
public:
  TranscriptPaths(LightGraph &,int maxPaths,int seqLen,const Model &);
  virtual ~TranscriptPaths();
  int numPaths() const;
  TranscriptPath *operator[](int);
  void computePosteriors();
  void computeLRs(double denom);
  void filter(double minScore);
  void sort(); // sorts descending by score
protected:
  LightGraph &G;
  int seqLen;
  const Model &model;
  Vector<TranscriptPath*> paths;
  void buildPaths(int maxPaths);
};

#endif

