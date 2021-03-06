/****************************************************************
 ACEplus_Vertex.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ACEplus_Vertex_H
#define INCL_ACEplus_Vertex_H
#include <iostream>
#include "LightVertex.H"
using namespace std;
using namespace BOOM;


class ACEplus_Vertex : public LightVertex {
public:
  ACEplus_Vertex(const String &substrate,SignalType,int begin,int end,
	      double score,Strand,int ID);
  void setRefScore(double);
  void setThreshold(double);
  void setSeq(const String &);
  double getRefScore() const;
  double getThreshold() const;
  const String &getSeq() const;
  void setRawScore(double);
  double getRawScore() const;
  void setDeNovo();
  bool isDeNovo() const;
protected:
  double refScore;  // Score of the corresponding signal in the reference
  double threshold; // Signal threshold
  double rawScore;  // probability before normalizing into LLR
  String seq;
  bool denovo;
};

#endif

