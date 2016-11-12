/****************************************************************
 Accumulator.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_Accumulator_H
#define INCL_Accumulator_H

#include "Propagator.H"

class Accumulator : public Propagator
{
public:
  inline Accumulator(int position);
};


inline Accumulator::Accumulator(int position) 
  : Propagator(position) 
{
  double *sc=scores;
  *sc=*(sc+1)=*(sc+2)=0;
}


#endif
