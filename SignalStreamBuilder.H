/****************************************************************
 SignalStreamBuilder.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalStreamBuilder_H
#define INCL_SignalStreamBuilder_H
#include <iostream>
#include "BOOM/Set.H"
#include "ReferenceAnnotation.H"
#include "VariantEvents.H"
#include "ConstraintIntervals.H"
#include "SignalStream.H"
#include "Isochore.H"
#include "SignalTypeProperties.H"
using namespace std;
using namespace BOOM;


class SignalStreamBuilder {
public:
  SignalStreamBuilder(const ReferenceAnnotation &, // input
		      const VariantEvents &,       // input
		      SignalStream &,              // output
		      ConstraintIntervals &,       // output
		      Set<Signal*> &newSignals,    // output
		      Isochore *,
		      int maxIntronScan,
		      int minExonLength,
		      int minIntronLen,
		      bool allowGains,
		      bool allowGainExonBrokenStop);
protected:
  const ReferenceAnnotation &refAnno;
  const VariantEvents &events;
  SignalStream &stream;
  ConstraintIntervals &constraints;
  Set<Signal*> &newSignals;
  Isochore *isochore;
  bool allowGains, allowGainExonBrokenStop;
  int maxIntronScan;
  int minExonLength;
  int minIntronLen;
  void build();
  void gain(const VariantEvent &);
  void loss(const VariantEvent &);
  void gainATG(int pos,SignalSensor &);
  void gainTAG(int pos,SignalSensor &);
  void gainGT(int pos,SignalSensor &);
  void gainAG(int pos,SignalSensor &);
  void lossATG(int pos,SignalSensor &);
  void lossTAG(int pos,SignalSensor &);
  void lossGT(int pos,SignalSensor &);
  void lossAG(int pos,SignalSensor &);
  void scan(int begin,int end,SignalSensor &);
};

#endif

