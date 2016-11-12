/****************************************************************
 MarkovChainCompiler.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_MarkovChainCompiler_H
#define INCL_MarkovChainCompiler_H

#include "MarkovChain.H"
#include "FastMarkovChain.H"

class MarkovChainCompiler
{
  static FastMarkovChain *compileForward(ContentSensor &);
  static FastMarkovChain *compileReverse(ContentSensor &);
public:
  static FastMarkovChain *compile(ContentSensor &MC_or_IMM_non_periodic);
  static BOOM::String stateToString(int state,FastMarkovChain &);
  static int stringToState(const BOOM::String &,FastMarkovChain &);
};


#endif
