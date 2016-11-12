/**************************************************************
 NthOrderStringIterator.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_NthOrderStringIterator_H
#define INCL_NthOrderStringIterator_H

#include "BOOM/Alphabet.H"
#include "BOOM/Sequence.H"
#include "BOOM/String.H"

// Generates strings of length N (not N+1!)

class NthOrderStringIterator
{
  Alphabet &alphabet;
  int alphabetSize, N;
  Sequence current;
  bool hasMore;

  void reset();
  void advance();
public:
  NthOrderStringIterator(int N,Alphabet &);
  BOOM::String getNextString();
  bool done();
};


#endif
