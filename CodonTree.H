/****************************************************************
 CodonTree.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_CodonTree_H
#define INCL_CodonTree_H
#include "BOOM/Array1D.H"
#include <iostream>
using namespace std;

class CodonTree
{
  BOOM::Array1D<double> array; // 125 elements (5^3 codons including 'N')

  void load(istream &);
  int codonToIndex(const char *);
  BOOM::String decode(int base5); // base-5 number -> codon
public:
  CodonTree(istream &);
  CodonTree(const CodonTree &);
  float scoreCodon(const char *,int pos);
  float scoreCodon(const char *);
  void revComplementSelf();
  CodonTree *clone() const;
};


#endif
