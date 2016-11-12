/****************************************************************
 GarbageCollector.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "GarbageCollector.H"
#include "Edge.H"
#include "Signal.H"
#include <iostream>

#undef USE_SMART_POINTERS

GarbageCollector::~GarbageCollector()
{
  // dtor

#ifdef EXPLICIT_GRAPHS
  BOOM::Set<SignalPtr>::iterator cur=reachableSignals.begin(), 
    end=reachableSignals.end();
  for(; cur!=end ; ++cur) delete *cur;
#endif
}



void GarbageCollector::addSignal(SignalPtr s)
{
#ifdef EXPLICIT_GRAPHS
  unreachableSignals.insert(s);
#endif
}



void GarbageCollector::markRight(SignalPtr origin)
{
#ifdef EXPLICIT_GRAPHS
  SignalPtr s=origin;
  if(!signalsReachableFromRight.isMember(origin)) return;
  signalsReachableFromRight.erase(origin);
  reachableSignals.insert(origin);

  BOOM::Set<Edge*> &edgesOut=s->getEdgesOut();
  BOOM::Set<Edge*>::iterator oCur=edgesOut.begin(), oEnd=edgesOut.end();
  for(; oCur!=oEnd ; ++oCur)
    {
      Edge *edge=*oCur;
      SignalPtr r=edge->getRight();
      markRight(r);
    }
#endif
}



void GarbageCollector::markLeft(SignalPtr origin)
{
#ifdef EXPLICIT_GRAPHS
  SignalPtr s=origin;
  if(!unreachableSignals.isMember(origin)) return;
  unreachableSignals.erase(origin);
  signalsReachableFromRight.insert(origin);

  BOOM::Set<Edge*> &edgesIn=s->getEdgesIn();
  BOOM::Set<Edge*>::iterator iCur=edgesIn.begin(), iEnd=edgesIn.end();
  for(; iCur!=iEnd ; ++iCur) {
      Edge *edge=*iCur;
      SignalPtr l=edge->getLeft();
      markLeft(l);
    }
#endif
}



void GarbageCollector::makeImmortal(SignalPtr signal)
{
  reachableSignals.insert(signal);
  unreachableSignals.erase(signal);
  signalsReachableFromRight.erase(signal);
}



void GarbageCollector::sweep()
{
#ifdef EXPLICIT_GRAPHS
  cerr<<"Collecting & deleting garbage..."<<endl;
  int R=reachableSignals.size();
  int U=unreachableSignals.size();
  int L=signalsReachableFromRight.size();
  int N=R+U+L;
  float percentDeleted=int((U+L)/float(N)*100+5/9.0);
  BOOM::Set<SignalPtr>::iterator cur=unreachableSignals.begin(), 
    end=unreachableSignals.end();
  for(; cur!=end ; ++cur) delete *cur;
  cur=signalsReachableFromRight.begin(), 
    end=signalsReachableFromRight.end();
  for(; cur!=end ; ++cur) delete *cur;
  unreachableSignals.clear();
  signalsReachableFromRight.clear();
  cerr<<percentDeleted<<"% of "<<N<<" vertices deleted, "<<R<<" remaining"<<endl;
#endif
}


void GarbageCollector::purge()
{
#ifdef EXPLICIT_GRAPHS
  BOOM::Set<SignalPtr>::iterator cur=unreachableSignals.begin(), 
    end=unreachableSignals.end();
  for(; cur!=end ; ++cur) delete *cur;
  cur=reachableSignals.begin(), end=reachableSignals.end();
  for(; cur!=end ; ++cur) delete *cur;
  cur=signalsReachableFromRight.begin(), end=signalsReachableFromRight.end();
  for(; cur!=end ; ++cur) delete *cur;
  unreachableSignals.clear();
  signalsReachableFromRight.clear();
  reachableSignals.clear();
#endif
}



BOOM::Set<SignalPtr>::iterator GarbageCollector::signalsBegin()
{
#ifdef EXPLICIT_GRAPHS
  return reachableSignals.begin();
#else
  throw "EXPLICIT_GRAPHS is not defined in GarbageColletor::signalsBegin()";
#endif
}



BOOM::Set<SignalPtr>::iterator GarbageCollector::signalsEnd()
{
#ifdef EXPLICIT_GRAPHS
  return reachableSignals.end();
#else
  throw "EXPLICIT_GRAPHS is not defined in GarbageColletor::signalsEnd()";
#endif
}



void GarbageCollector::drop(SignalPtr s)
{
#ifdef EXPLICIT_GRAPHS
  if(reachableSignals.isMember(s)) 
    reachableSignals.remove(s);
  else if(unreachableSignals.isMember(s)) 
    unreachableSignals.remove(s);
#endif
}
