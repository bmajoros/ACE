/****************************************************************
 RnaJunction.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_RnaJunction_H
#define INCL_RnaJunction_H
#include <math.h>
#include <iostream>
#include "BOOM/File.H"
#include "BOOM/Comparator.H"
#include "BOOM/Strand.H"
using namespace std;
using namespace BOOM;

class JunctionComparator;

class RnaJunction {
public:
  RnaJunction();
  RnaJunction(int b,int e,float d,char s) : begin(b), end(e), depth(d), 
					  strand(s) {}
  bool read(File &);
  bool write(File &);
  int getBegin();
  int getEnd();
  float getDepth();
  void shift(int add) { begin+=add; end+=add; }
  Strand getStrand();
  static int getSize() { return 3*sizeof(int)+sizeof(char); }
  void printOn(ostream &) const;
  void logify() { depth=log(depth); }
protected:
  int begin, end;
  float depth;
  Strand strand;
  friend class JunctionComparator;
  friend class JunctionEndComparator;
};

ostream &operator<<(ostream &,const RnaJunction &);

class JunctionComparator : public Comparator<RnaJunction> {
  bool equal(RnaJunction &a,RnaJunction &b)   {return a.begin==b.begin;}
  bool greater(RnaJunction &a,RnaJunction &b) {return a.begin>b.begin;}
  bool less(RnaJunction &a,RnaJunction &b)    {return a.begin<b.begin;}
};

class JunctionEndComparator : public Comparator<RnaJunction> {
  bool equal(RnaJunction &a,RnaJunction &b)   {return a.end==b.end;}
  bool greater(RnaJunction &a,RnaJunction &b) {return a.end>b.end;}
  bool less(RnaJunction &a,RnaJunction &b)    {return a.end<b.end;}
};

#endif
