/****************************************************************
 ContentRegions.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ContentRegions_H
#define INCL_ContentRegions_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/GffTranscript.H"
#include "ContentRegion.H"
using namespace std;
using namespace BOOM;

class ContentRegions {
public:
  ContentRegions(const GffTranscript &,int sequenceLength);
  const Vector<ContentRegion> &asVector();
  bool findJunction(int pos,const ContentRegion *&preceding,
		    const ContentRegion *&following) const;
  ContentRegion *regionOverlapping(int pos) const;
  void printOn(ostream &) const;
protected:
  Vector<ContentRegion> regions;
  void init(const GffTranscript &,int seqLen);
  void sort();
};

ostream &operator<<(ostream &,const ContentRegions &);

#endif

