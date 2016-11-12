/****************************************************************
 SignalLabelingProfile.H
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalLabelingProfile_H
#define INCL_SignalLabelingProfile_H
#include <iostream>
#include "BOOM/Array2D.H"
#include "SignalSensor.H"
#include "Labeling.H"
using namespace std;
using namespace BOOM;

/* This class represents the individual nucleotide labelings that
   would occur if a signal was to occur in a particular phase.  It
   provides precomputed labelings for all positions in the PWM (signal
   context window), for constant-time lookup when computing priors
   under a given reference labeling.
 */
class SignalLabelingProfile {
public:
  SignalLabelingProfile(SignalSensor &);
  SignalLabelingProfile() {}
  GeneModelLabel getLabel(int signalPhase,int windowPos);
  void printOn(ostream &);
  int getLength() const;
protected:
  Array2D<GeneModelLabel> M; // indexed as: [phase][position]
  void init(SignalSensor &);
  void initATG(int offset,int len);
  void initTAG(int offset,int len);
  void initGT(int offset,int len);
  void initAG(int offset,int len);
  void initNegATG(int offset,int len);
  void initNegTAG(int offset,int len);
  void initNegGT(int offset,int len);
  void initNegAG(int offset,int len);
};
ostream &operator<<(ostream &,const SignalLabelingProfile &);

#endif

