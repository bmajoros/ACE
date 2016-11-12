/****************************************************************
 SignalLabelingProfile.C
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SignalLabelingProfile.H"
using namespace std;
using namespace BOOM;

SignalLabelingProfile::SignalLabelingProfile(SignalSensor &ss)
{
  init(ss);
}



GeneModelLabel SignalLabelingProfile::getLabel(int signalPhase,int windowPos)
{
  return M[signalPhase][windowPos];
}



void SignalLabelingProfile::init(SignalSensor &ss)
{
  const int L=ss.getContextWindowLength();
  SignalType t=ss.getSignalType();
  int offset=ss.getConsensusOffset();
  M.resize(3,L);
  switch(t)
    {
    case ATG:       initATG(offset,L); break;
    case TAG:       initTAG(offset,L); break;
    case GT:        initGT(offset,L); break;
    case AG:        initAG(offset,L); break;
    case NEG_ATG:   initNegATG(offset,L); break;
    case NEG_TAG:   initNegTAG(offset,L); break;
    case NEG_GT:    initNegGT(offset,L); break;
    case NEG_AG:    initNegAG(offset,L); break;
    default:        M.setAllTo(LABEL_NONE); break;
    }
}



void SignalLabelingProfile::initATG(int offset,int len)
{
  //----ATG----   => signal phase is 0
  //    0120120
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset, phase=0 ; i<len ; ++i, phase=(phase+1)%3)
    M[0][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initTAG(int offset,int len)
{
  //----TAG----   => signal phase is 0
  //2012012
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset+2, phase=2 ; i>=0 ; --i, phase=posmod(phase-1))
    M[0][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initGT(int offset,int len)
{
  //----GT----    => signal phase is 1
  //0120
  M.setAllTo(LABEL_INTRON);
  for(int signalPhase=0 ; signalPhase<3 ; ++signalPhase)
    for(int i=offset-1, phase=posmod(signalPhase-1) ; i>=0 ; 
	--i, phase=posmod(phase-1))
      M[signalPhase][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initAG(int offset,int len)
{
  //----AG----    => signal phase is 2
  //      2012
  M.setAllTo(LABEL_INTRON);
  for(int signalPhase=0 ; signalPhase<3 ; ++signalPhase)
    for(int i=offset+2, phase=signalPhase ; i<len ; ++i, phase=(phase+1)%3)
      M[signalPhase][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initNegATG(int offset,int len)
{
  //----CAT----   => signal phase is 2 (reverse-strand start codon)
  //0210210
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset+2, phase=0 ; i>=0 ; --i, phase=(phase+1)%3)
    M[0][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initNegTAG(int offset,int len)
{
  //----CTA----   => signal phase is 2 (reverse-strand stop codon)
  //    2102102
  M.setAllTo(LABEL_INTERGENIC);
  for(int i=offset, phase=2 ; i<len ; ++i, phase=posmod(phase-1))
    M[0][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initNegGT(int offset,int len)
{
  //----AC----    => signal phase is 0 (reverse-strand donor)
  //      0210
  M.setAllTo(LABEL_INTRON);
  for(int signalPhase=0 ; signalPhase<3 ; ++signalPhase)
    for(int i=offset+2, phase=signalPhase ; i<len ; ++i, phase=posmod(phase-1))
	M[signalPhase][i]=getExonLabel(phase);
}



void SignalLabelingProfile::initNegAG(int offset,int len)
{
  //----CT----    => signal phase is 1 (reverse-strand acceptor)
  //2102
  M.setAllTo(LABEL_INTRON);
  for(int signalPhase=0 ; signalPhase<3 ; ++signalPhase)
    for(int i=offset-1, phase=signalPhase ; i>=0 ; --i, phase=posmod(phase-1))
	M[signalPhase][i]=getExonLabel(phase);
}



int SignalLabelingProfile::getLength() const
{
  return M.getSecondDim();
}



void SignalLabelingProfile::printOn(ostream &os)
{
  int L=getLength();
  for(int phase=0 ; phase<3 ; ++phase) {
    for(int pos=0 ; pos<L ; ++pos) os<<M[phase][pos]<<" ";
    os<<endl;
  }
}



ostream &operator<<(ostream &os,const SignalLabelingProfile &p)
{
  p.printOn(os);
  return os;
}



