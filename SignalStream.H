/****************************************************************
 SignalStream.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalStream_H
#define INCL_SignalStream_H
#include <iostream>
#include "BOOM/Vector.H"
#include "Signal.H"
using namespace std;
using namespace BOOM;


class SignalStream {
public:
  SignalStream();
  virtual ~SignalStream();
  void add(Signal *);
  void sort();
  void deduplicate();
  const Vector<Signal*> &peek() { return signals; }
  void reset();
  Signal *detect(int position);
  void printOn(ostream &) const;
protected:
  Vector<Signal*> signals;
  int currentIndex;
};

ostream &operator<<(ostream &,const SignalStream &);

#endif

