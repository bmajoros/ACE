/****************************************************************
 SignalSensors.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalSensors_H
#define INCL_SignalSensors_H
#include <iostream>
#include "BOOM/Set.H"
#include "BOOM/String.H"
#include "SignalSensor.H"
using namespace std;
using namespace BOOM;


struct SignalSensors {
  SignalSensor *startCodonSensor;
  SignalSensor *shortStartSensor; // just the ATG
  SignalSensor *stopCodonSensor;
  SignalSensor *donorSensor;
  SignalSensor *acceptorSensor;
  Set<String> donorConsensuses;
  Set<String> acceptorConsensuses;
  Set<String> startCodons;
  Set<String> stopCodons;
  SignalSensors();
  void setConsensuses();
  SignalSensor *findSensor(SignalType);
private:
  void setConsensuses(SignalSensor *,const Set<String> &consensuses);
};

#endif

