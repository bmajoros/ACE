/****************************************************************
 SignalSensors.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SignalSensors.H"
using namespace std;
using namespace BOOM;

SignalSensors::SignalSensors()
  : startCodonSensor(NULL), stopCodonSensor(NULL), shortStartSensor(NULL),
    donorSensor(NULL), acceptorSensor(NULL)
{
  //ctor
}



void SignalSensors::setConsensuses()
{
  setConsensuses(startCodonSensor,startCodons);
  setConsensuses(shortStartSensor,startCodons);
  setConsensuses(stopCodonSensor,stopCodons);
  setConsensuses(donorSensor,donorConsensuses);
  setConsensuses(acceptorSensor,acceptorConsensuses);
}



void SignalSensors::setConsensuses(SignalSensor *sensor,
				   const Set<String> &consensuses)
{
  StringMap<char> &known=sensor->getConsensuses();
  for(Set<String>::const_iterator cur=consensuses.begin(), end=
	consensuses.end() ; cur!=end ; ++cur) {
    const String &s=*cur;
    if(!known.isDefined(s)) sensor->addConsensus(s);
  }
}



SignalSensor *SignalSensors::findSensor(SignalType type)
{
  switch(type) {
  case ATG: return startCodonSensor;
  case TAG: return stopCodonSensor;
  case GT:  return donorSensor;
  case AG:  return acceptorSensor;
  case TSS: return NULL;
  case TES: return NULL;
  case LEFT_TERMINUS:   return NULL;
  case RIGHT_TERMINUS:  return NULL;
  }
  cout<<"No sensor found for signal type "<<signalTypeToString(type)<<endl;
  INTERNAL_ERROR;
}



