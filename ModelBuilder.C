/**************************************************************
 ModelBuilder.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include "ModelBuilder.H"
#include <iostream>
#include "MarkovChain.H"
#include "WMM.H"
#include "WAM.H"
#include "WWAM.H"
#include "ThreePeriodicMarkovChain.H"
#include "IMM.H"
#include "ThreePeriodicIMM.H"
#include "BranchAcceptor.H"


ModelBuilder::ModelBuilder(GarbageCollector &gc,Alphabet &alphabet,
			   int minSampleSize,int order,int windowSize)
  : alphabet(alphabet), 
    minSampleSize(minSampleSize), 
    order(order),
    windowSize(windowSize),
    gc(gc),
    branchContextLength(19) // default for GENSCAN is 19
{
  // ctor
}



ContentSensor *ModelBuilder::buildContentSensor(ModelType modelType,
					BOOM::Vector<TrainingSequence*> &seqs,
					ContentType contentType)
{
  switch(modelType)
    {
    case MARKOV_CHAIN:
      return new MarkovChain(seqs,order,minSampleSize,
			     MarkovChain::NO_PHASE,contentType);
    case THREE_PERIODIC:
      return new ThreePeriodicMarkovChain(seqs,order,minSampleSize,
					  contentType);
    case IMM_MODEL:
      return new IMM(seqs,order,minSampleSize,
		     MarkovChain::NO_PHASE,contentType);
    case IMM_3P:
      return new ThreePeriodicIMM(seqs,order,minSampleSize,contentType);
    }
}



SignalSensor *ModelBuilder::buildSignalSensor(ModelType modelType,
				      BOOM::Vector<TrainingSequence*> &
					      sequences,
					      SignalType signalType,
					      int consensusOffset,
					      int consensusLength)
{
  switch(modelType)
    {
    case WMM_MODEL:
      return new WMM(gc,sequences,signalType,consensusOffset,
		     consensusLength);
    case WAM_MODEL:
      return new WAM(gc,sequences,order,minSampleSize,signalType,
		     consensusOffset,consensusLength);
    case WWAM_MODEL:
      return new WWAM(gc,sequences,order,minSampleSize,signalType,
		      consensusOffset,consensusLength);

    case BRANCH_ACCEPTOR:
      return new BranchAcceptor(gc,sequences,branchOrder,order,
				branchContextLength,minSampleSize,
				consensusOffset,AG);

      /*
    case MARKOV_CHAIN:
      return new NthOrderMarkovChain(sequences,order,alphabet,minSampleSize);
    case IMM_MODEL:
      return new IMM(sequences,order,alphabet,minSampleSize);
    case CODON_BIAS:
      return new CodingPotential(sequences,alphabet);
      */

    default:
      throw "Error in ModelBuilder::buildModel()";
    }
}



void ModelBuilder::changeOrder(int order)
{
  this->order=order;
}



void ModelBuilder::changeWindowSize(int windowSize)
{
  this->windowSize=windowSize;
}



void ModelBuilder::changeMinSampleSize(int minSampleSize)
{
  this->minSampleSize=minSampleSize;
}



void ModelBuilder::changeBranchContextLength(int l)
{
  branchContextLength=l;
}



void ModelBuilder::changeBranchOrder(int o)
{
  branchOrder=o;
}


