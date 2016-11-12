/****************************************************************
 RnaJunctions.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_RnaJunctions_H
#define INCL_RnaJunctions_H
#include <iostream>
#include "BOOM/File.H"
#include "BOOM/Vector.H"
#include "RnaJunction.H"
using namespace std;
using namespace BOOM;

// ### NOTE: WE CURRENTLY IGNORE STRAND INFORMATION!!!

// Precondition: the junction file is already sorted by begin coordinate

class RnaJunctions {
public:
  RnaJunctions(const String &filename="");
  bool load(const String &filename);
  int getNumJunctions() const;
  const RnaJunction &operator[](int) const;
  float getDepth(int begin,int end) const; // zero-based, half-open
  float getSpliceOutDepth(int pos) const; // leftmost intron pos
  float getSpliceInDepth(int pos) const; // rightmost intron pos + 1
  float getMaxDepth() const;
  void logify();
protected:
  Vector<RnaJunction> junctions;
  Vector<int> sortedByEnd;
  void sortEnds();
};

#endif
