/****************************************************************
 EdgeFactory.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "EdgeFactory.H"
#include <iostream>
#include "Signal.H"
#include "SignalTypeProperties.H"


/****************************************************************
                      EdgeFactory methods
 ****************************************************************/

PhasedEdge *EdgeFactory::newPhasedEdge(double scorePhase0,double scorePhase1,
				       double scorePhase2,SignalPtr left,
				       SignalPtr right)
{
  return new PhasedEdge(scorePhase0,scorePhase1,scorePhase2,left,right);
}



NonPhasedEdge *EdgeFactory::newNonPhasedEdge(double score,SignalPtr left,
					     SignalPtr right)
{
  return new NonPhasedEdge(score,left,right);
}



/****************************************************************
                  FilteredEdgeFactory methods
 ****************************************************************/

FilteredEdgeFactory::FilteredEdgeFactory(EvidenceFilter *f)
  : filter(f)
{
  // ctor
}



FilteredEdgeFactory::~FilteredEdgeFactory() {
  delete filter; 
}



PhasedEdge *FilteredEdgeFactory::newPhasedEdge(double scorePhase0,
					       double scorePhase1,
					       double scorePhase2,
					       SignalPtr left,
					       SignalPtr right)
{
  SignalType leftType=left->getSignalType(), rightType=right->getSignalType();
  ContentType T=
    SignalTypeProperties::global.getContentType(leftType,rightType);
  int leftPos=left->getConsensusPosition();
  int rightPos=right->getConsensusPosition();
  if(isIntron(T)) {
    //cout<<"******* "<<leftPos<<" "<<rightPos+2<<endl;
    if(!filter->intronSupported(leftPos,rightPos+2)) {
      //cout<<"\tREJECTING INTRON "<<leftPos<<"-"<<rightPos+2<<endl;
      return NULL;
    }
    //else cout<<"supporting intron "<<leftPos<<"-"<<rightPos+2<<endl;
  }
  else if(isCoding(T)) {
    if(endsIntron(leftType)) leftPos+=2;
    if(beginsIntron(rightType)) rightPos-=3;

    //cout<<"CHECKING PILEUP FOR INTERVAL "<<leftPos<<"-"<<rightPos+3<<endl;


    //if(rightType==TAG) cout<<"checking exon "<<leftPos<<" to "<<rightPos+3<<endl;
    //if(filter->exonSupported(leftPos,rightPos+3)) return NULL;


  }
  else INTERNAL_ERROR; // intergenic -- should not happen
  //if(T==INTRON) cout<<"INTRON EDGE "<<leftPos<<"-"<<rightPos<<endl;
  return new PhasedEdge(scorePhase0,scorePhase1,scorePhase2,left,right);
}



void FilteredEdgeFactory::setEvidence(EvidenceFilter *f)
{
  delete filter;
  filter=f;
}


