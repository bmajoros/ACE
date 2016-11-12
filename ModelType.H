/*
 ModelType.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#ifndef INCL_ModelType_H
#define INCL_ModelType_H

#include "BOOM/String.H"

enum ModelType
  {
    WMM_MODEL,        // simple weight matrix (assumes independence)
    WAM_MODEL,        // Weight Array Matrix (WMM w/ position dependencies)
    WWAM_MODEL,       // Windowed WAM (averages cond. probs. over a window)
    THREE_PERIODIC,   // 3-periodic (inhomogeneous) Markov Chain
    MARKOV_CHAIN,     // homogeneous Markov Chain
    FAST_MC,          // optimized markov chain
    FAST_3PMC,        // optimized 3-periodic markov chain
    IMM_MODEL,        // Interpolated Markov Model (ave. over Markov orders)
    CODON_BIAS,       // codon bias
    MDD,              // Maximal Dependence Decomposition tree
    SIGNAL_MODEL,     // An MDD tree with left and right context models
    HMM_MODEL,        // Simple HMM; each state emits one nucleotide
    NONSTATIONARY_MC, // Nonstationary Markov Chain (NSMC)
    IMM_3P,           // 3-period interpolated Markov chain
    BRANCH_ACCEPTOR,  // branch point WWAM + acceptor site WAM, like GENSCAN
  };

BOOM::String toString(ModelType);
ModelType stringToModelType(const BOOM::String &);

class ModelTypeInitializer
{
  static ModelTypeInitializer initializer;
  ModelTypeInitializer();
};


#endif
