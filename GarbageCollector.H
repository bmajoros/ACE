/****************************************************************
 GarbageCollector.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_GarbageCollector_H
#define INCL_GarbageCollector_H
#include "genezilla.H"
#include "Signal.H"

#undef USE_SMART_POINTERS

/****************************************************************
  A NOTE ON GARBAGE COLLECTION VS. SMART-POINTERS

  Deleting unused objects can be done by either (1) directly deleting
  them when it is clear they are no longer needed, using the delete
  operator, (2) periodically collecting objects that are no longer in
  use and deleting them in batches, or (3) using a "smart-pointer"
  mechanism to cause objects to auto-delete when their reference count
  drops to zero.

  When GeneZilla is run in USE_EXPLICIT_GRAPHS mode, we use garbage
  collection (GC) because the GC here actually eliminates not only
  those Signal objects no longer reachable, but also those that are
  reachable from one end of the graph but not the other -- this because
  the graph represents complete parses of the sequence.  A signal
  reachable from only one end of the graph cannot participate in a 
  full parse of the sequence.

  When GeneZilla is run in the non-graph mode, a "trellis" is built
  instead of the full graph.  The trellis consists of only a single
  pointer at each signal (per phase) to the optimal predecessor.
  In this case, GC is not needed, and indeed we can reduce the
  instantaneous memory consumption by using smart pointers rather than
  waiting for memory to become potentially very cluttered before the
  GC is run.  

  Thus, we observe this rule:

     if using graphs -> then employ garbage collection
     otherwise -> employ smart pointers

  Note that these two mechanisms must NOT be used together!
 ****************************************************************/


class GarbageCollector
{
protected:
  BOOM::Set<SignalPtr> signalsReachableFromRight; // reachable from right
  BOOM::Set<SignalPtr> reachableSignals; // reachable from left & right
  BOOM::Set<SignalPtr> unreachableSignals; // not reachable from both left & right
public:
  virtual ~GarbageCollector();
  virtual void addSignal(SignalPtr);
  void makeImmortal(SignalPtr);
  void markLeft(SignalPtr rightTerminus); // do this first! (restores order)
  void markRight(SignalPtr leftTerminus); // do this second! (reverses order)
  void sweep(); // delete unmarked signals
  BOOM::Set<SignalPtr>::iterator signalsBegin();
  BOOM::Set<SignalPtr>::iterator signalsEnd();
  BOOM::Set<SignalPtr>::iterator unreachableSignalsBegin() 
    {return unreachableSignals.begin();}
  BOOM::Set<SignalPtr>::iterator unreachableSignalsEnd() 
    {return unreachableSignals.end();}
  void drop(SignalPtr s);
  void purge(); // delete everything
};


// A GarbageCollector that ignores garbage given to it -- useful when
// allocating temporary Signals that require a GC parameter but aren't
// safe to give to the real GC.
class GarbageIgnorer : public GarbageCollector
{
public:
  virtual void addSignal(SignalPtr s) {reachableSignals.insert(s);}
};

#endif
