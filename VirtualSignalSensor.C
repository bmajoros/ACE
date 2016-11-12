/****************************************************************
 VirtualSignalSensor.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VirtualSignalSensor.H"
using namespace std;
using namespace BOOM;


VirtualSignalSensor::VirtualSignalSensor(SignalType t,GarbageCollector &gc)
  : SignalSensor(gc)
{
  setStrand(FORWARD_STRAND);
  setSignalType(t);
  setSizes(0,0,0);
}



bool VirtualSignalSensor::consensusOccursAt(const BOOM::String &,int index)
{
  return true;
}



SignalPtr VirtualSignalSensor::detect(const Sequence &,const BOOM::String &,
				      int contextWindowPosition)
{
  return new Signal(contextWindowPosition,0.0,*this,getGC(),getSignalType());
}



SignalPtr VirtualSignalSensor::detectWithNoCutoff(const Sequence &,
			     const BOOM::String &,
			     int contextWindowPosition)
{
  return new Signal(contextWindowPosition,0.0,*this,getGC(),getSignalType());
}



SignalSensor *VirtualSignalSensor::reverseComplement()
{
  return this;
}



double VirtualSignalSensor::getLogP(const Sequence &,const BOOM::String &,int begin)
{
  return 0;
}



bool VirtualSignalSensor::save(const BOOM::String &filename)
{
  return false;
}



bool VirtualSignalSensor::save(ostream &os)
{
  return false;
}



void VirtualSignalSensor::useLogOdds(SignalSensor &nullModel)
{
}



void VirtualSignalSensor::useLogOdds_anonymous(ContentSensor &nullModel)
{
}


