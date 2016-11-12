/****************************************************************
 Transitions.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_Transitions_H
#define INCL_Transitions_H

#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Array2D.H"
#include "SignalType.H"

class Transitions
{
  BOOM::Array2D<double> matrix;
  float optimism, intronOptimism;
  void load(istream &);
  void normalize();
  void logify();
public:
  Transitions(int numSignalTypes,istream &,float optimism,
	      float intronOptimism);
  double getLogP(SignalType from,SignalType to);
};


#endif
