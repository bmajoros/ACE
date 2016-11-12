/****************************************************************
 SignalTypeProperties.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalTypeProperties_H
#define INCL_SignalTypeProperties_H

#include "BOOM/Set.H"
#include "BOOM/Map.H"
#include "genezilla.H"
#include "BOOM/Strand.H"
#include "SignalType.H"
#include "ContentType.H"

class SignalTypeProperties
{
  struct SignalTypeInfo {
    BOOM::Set<int> allowablePhases;
    BOOM::Set<ContentType> belongsTo; // queues it can belong to
    BOOM::Set<ContentType> linksBackTo; // queues it can link back to
    bool codingConsensus;
    BOOM::Map<ContentType,int> propagatorIndexMap;
    Strand strand;
    bool enabled;
    SignalTypeInfo() : codingConsensus(false), strand(FORWARD_STRAND),
		       enabled(false) {}
  };
  BOOM::Map<SignalType,SignalTypeInfo> info;
public:
  static SignalTypeProperties global;
  SignalTypeProperties();
  bool hasBeenLoaded() {return !info.isEmpty();}

  // Storing information about SignalTypes:
  SignalTypeProperties &allowPhase(SignalType,int phase);
  SignalTypeProperties &belongTo(SignalType,ContentType queueType);
  SignalTypeProperties &linkTo(SignalType,ContentType queueType);
  void setStrand(SignalType,Strand);
  void setConsensusCoding(SignalType,bool);
  void enable(SignalType);
  void disable(SignalType);

  // Querying information about SignalTypes:
  inline BOOM::Set<int> &getAllowablePhases(SignalType);
  inline BOOM::Set<ContentType> &belongsInWhichQueues(SignalType);
  inline BOOM::Set<ContentType> &linksBackToWhichQueues(SignalType);
  inline int whichPropagator(SignalType,ContentType);
  inline bool consensusBasesAreCoding(SignalType);
  inline Strand getStrand(SignalType);
  inline ContentType getContentType(SignalType from,SignalType to);
  inline bool isEnabled(SignalType);
};



inline BOOM::Set<int> &SignalTypeProperties::getAllowablePhases(SignalType t)
{
  if(!info.isDefined(t)) throw "Undefined signal type in SignalTypeProperties::getAllowablePhases(): "+signalTypeToString(t);
  return info[t].allowablePhases;
}



inline BOOM::Set<ContentType> &SignalTypeProperties::belongsInWhichQueues(
   SignalType t)
{
  if(!info.isDefined(t)) throw "Undefined signal type in SignalTypeProperties::belongsInWhichQueus(): "+signalTypeToString(t);
  return info[t].belongsTo;
}



inline BOOM::Set<ContentType> &SignalTypeProperties::linksBackToWhichQueues(
   SignalType t)
{
  if(!info.isDefined(t)) throw "Undefined signal type in SignalTypeProperties::linksBackToWhichQueues(): "+signalTypeToString(t);
  return info[t].linksBackTo;
}



inline int SignalTypeProperties::whichPropagator(SignalType t,ContentType c)
{
  if(!info.isDefined(t)) throw "Undefined signal type in SignalTypeProperties::whichPropagator(): "+signalTypeToString(t);
  return info[t].propagatorIndexMap[c];
}



inline bool SignalTypeProperties::consensusBasesAreCoding(SignalType t)
{
  if(!info.isDefined(t)) throw "Undefined signal type in SignalTypeProperties::consensusBasesAreCoding(): "+signalTypeToString(t);
  return info[t].codingConsensus;
}



inline Strand SignalTypeProperties::getStrand(SignalType t)
{
  if(!info.isDefined(t)) throw "Undefined signal type in SignalTypeProperties::getStrand(): "+signalTypeToString(t);
  return info[t].strand;
}



inline bool SignalTypeProperties::isEnabled(SignalType t)
{
  if(!info.isDefined(t)) throw "Undefined signal type in SignalTypeProperties::isEnabled(): "+signalTypeToString(t);
  return info[t].enabled;
}



inline ContentType SignalTypeProperties::getContentType(SignalType from,
							SignalType to)
{
  BOOM::Set<ContentType> lType=belongsInWhichQueues(from);
  BOOM::Set<ContentType> rType=linksBackToWhichQueues(to);
  BOOM::Set<ContentType> cType;
  lType.intersectWith(rType,cType);
  return *cType.begin();
}



#endif
