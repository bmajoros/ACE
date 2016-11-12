/****************************************************************
 Edge.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "Edge.H"
#include <iostream>
#include "Signal.H"
#include "SignalTypeProperties.H"
#include "Propagator.H"
#include "genezilla.H"


Edge::Edge(SignalPtr left,SignalPtr right)
  : left(left),
    right(right)
{
}



Edge::~Edge()
{
}



SignalPtr Edge::getLeft()
{
  return left;
}



SignalPtr Edge::getRight()
{
  return right;
}



int Edge::getEdgeBegin()
{
  return left->getContextWindowPosition()+left->getContextWindowLength();
}



int Edge::getEdgeEnd()
{
  return right->getContextWindowPosition();
}



int Edge::getFeatureBegin()
{
  int consensusPos=left->getConsensusPosition();
  switch(left->getSignalType())
    {
    case LEFT_TERMINUS:  return 0;
    case RIGHT_TERMINUS: return left->getContextWindowPosition();
    case ATG:       return consensusPos;     // exon
    case TAG:       return consensusPos+3;   // UTR or intergenic
    case UTR5GT:
    case UTR3GT:
    case GT:        return consensusPos;     // intron
    case UTR5AG:
    case UTR3AG:
    case AG:        return consensusPos+2;   // exon
    case TSS:       return consensusPos;     // UTR
    case TES:       return consensusPos+left->getConsensusLength(); // intergenic
    case NEG_ATG:   return consensusPos+3;   // UTR or intergenic
    case NEG_TAG:   return consensusPos;     // exon
    case NEG_GT:    return consensusPos+2;   // exon
    case NEG_AG:    return consensusPos;     // intron
    case NEG_TSS:   return consensusPos+left->getConsensusLength(); // intergenic
    case NEG_TES:   return consensusPos;     // UTR
    default: INTERNAL_ERROR;
    }
}



int Edge::getFeatureEnd()
{
  int consensusPos=right->getConsensusPosition();
  switch(right->getSignalType())
    {
    case LEFT_TERMINUS:  return 0;
    case RIGHT_TERMINUS: return right->getContextWindowPosition();
    case ATG:       return consensusPos;     // UTR or intergenic
    case TAG:       return consensusPos+3;   // exon
    case UTR5GT:
    case UTR3GT:
    case GT:        return consensusPos;     // exon
    case UTR5AG:
    case UTR3AG:
    case AG:        return consensusPos+2;   // intron
    case TSS:       return consensusPos;     // intergenic
    case TES:       return consensusPos+right->getConsensusLength(); // UTR
    case NEG_ATG:   return consensusPos+3;   // exon
    case NEG_TAG:   return consensusPos;     // UTR or intergenic
    case NEG_GT:    return consensusPos+2;   // intron
    case NEG_AG:    return consensusPos;     // exon
    case NEG_TSS:   return consensusPos+right->getConsensusLength(); // UTR
    case NEG_TES:   return consensusPos;     // intergenic
    default: INTERNAL_ERROR;
    }
}



ContentType Edge::getContentType() const
{
  return 
    SignalTypeProperties::global.getContentType(
      left->getSignalType(),
      right->getSignalType());
}



Strand Edge::getStrand()
{
  return ::getStrand(getContentType());
}



bool Edge::isCoding()
{
  return ::isCoding(getContentType());
}



bool Edge::isIntron()
{
  return ::isIntron(getContentType());
}



bool Edge::isIntergenic()
{
  return ::isIntergenic(getContentType());
}



bool Edge::isUTR()
{
  return ::isUTR(getContentType());
}



int Edge::propagateForward(int phase)
{
  if(isIntergenic()) return (right->getStrand()==FORWARD_STRAND ? 0 : 2);
  if(!isCoding()) return phase;
  int length=getFeatureEnd()-getFeatureBegin();
  switch(getStrand())
    {
    case FORWARD_STRAND:
      return (phase+length)%3;
    case REVERSE_STRAND:
      return posmod(phase-length);
    }
}



int Edge::propagateBackward(int phase)
{
  if(isIntergenic()) return (left->getStrand()==FORWARD_STRAND ? 0 : 2);
  if(!isCoding()) return phase;
  int length=getFeatureEnd()-getFeatureBegin();
  switch(getStrand())
    {
    case FORWARD_STRAND:
      return posmod(phase-length);
    case REVERSE_STRAND:
      return (phase+length)%3;
    }
}



void Edge::printOn(ostream &os) const
{
  left->printOn(os);os<<" -> ";right->printOn(os);
  os<<" ("<<getContentType()<<")";
}



ostream &operator<<(ostream &os,const Edge &edge)
{
  edge.printOn(os);
  return os;
}



//=================================================================
//                        PhasedEdge methods
//=================================================================

PhasedEdge::PhasedEdge(double scorePhase0,double scorePhase1,
		       double scorePhase2,
		       SignalPtr left,SignalPtr right)
  : Edge(left,right)
{
  scores[0]=scorePhase0;
  scores[1]=scorePhase1;
  scores[2]=scorePhase2;
}



double PhasedEdge::getEdgeScore(int phase)
{
  return scores[phase];
}



void PhasedEdge::setEdgeScore(int phase,double s)
{
  scores[phase]=s;
}



bool PhasedEdge::isPhased() const
{
  return true;
}



double PhasedEdge::getInductiveScore(int phase)
{
  return scores[phase]+left->posteriorInductiveScore(phase);
}



void PhasedEdge::printOn(ostream &os) const
{
  os<<"[";
  for(int i=0 ; i<3 ; ++i)
    if(!isinf(scores[i])) {
      os<<scores[i];
      if(i<2 && !isinf(scores[i+1])) os<<',';
    }
  os<<"] ";
  Edge::printOn(os);
}



//=================================================================
//                      NonPhasedEdge methods
//=================================================================
NonPhasedEdge::NonPhasedEdge(double score,SignalPtr left,SignalPtr right)
  : score(score),
    Edge(left,right)
{
}



double NonPhasedEdge::getEdgeScore()
{
  return score;
}



void NonPhasedEdge::setEdgeScore(double s)
{
  score=s;
}



bool NonPhasedEdge::isPhased() const
{
  return false;
}



double NonPhasedEdge::getInductiveScore()
{
  // Invariant: contentType is INTERGENIC or *_UTR

  Strand strand=left->getStrand();
  int defaultPhase=(strand==FORWARD_STRAND ? 0 : 2);
  return score+left->posteriorInductiveScore(defaultPhase);
}



void NonPhasedEdge::printOn(ostream &os) const
{
  //os<<"[0] ";
  os<<"["<<score<<"] ";
  Edge::printOn(os);
}


