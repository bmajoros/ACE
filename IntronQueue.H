/****************************************************************
 IntronQueue.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_IntronQueue_H
#define INCL_IntronQueue_H
#include "BOOM/Map.H"
#include "BOOM/Set.H"
#include "SignalQueue.H"


class IntronQueueIterator : public BOOM::Iterator<SignalPtr>
{
  FixedSizePriorityList<SignalPtr >::iterator iterators[3];
  FixedSizePriorityList<SignalPtr >::iterator ends[3];
  int index;

  void increment();
public:
  IntronQueueIterator() {}
  IntronQueueIterator(FixedSizePriorityList<SignalPtr > *lists[3],
		      bool pastTheEnd=false);
  virtual void operator++();
  virtual void operator++(int);
  virtual SignalPtr &operator*();
  virtual bool operator==(const BOOM::Iterator<SignalPtr > &) const;
  virtual BOOM::Iterator<SignalPtr > &clone();
  FixedSizePriorityList<SignalPtr >::iterator getNativeIterator() const;
  int getIndex() const;
};



/*************************************************************************
  This is a special type of SignalQueue in which there can be multiple
  (overlapping) signals in the holding queue, but there can be at most
  three signals in the main queue -- a best signal for each phase.  When 
  a signal is ready to graduate from the holding queue to the main queue 
  we see (in each phase) if it is better than one in the main queue (i.e., 
  has a better inductive score in that phase), and if so, move it there; 
  otherwise we discard it.
**************************************************************************/

class IntronQueue : public SignalQueue
{
  class MembershipCounter // counts # "main" queues a signal is in (0-3)
  {
    BOOM::Map<SignalPtr ,int> counter;
  public: 
    void increment(SignalPtr s) 
      {if(counter.isDefined(s)) ++counter[s]; else counter[s]=1;}
    void decrement(SignalPtr s) 
      {--counter[s];if(!counter[s]) counter.undefine(s);}
    BOOM::Map<SignalPtr ,int>::iterator begin() {return counter.begin();}
    BOOM::Map<SignalPtr ,int>::iterator end() {return counter.end();}
    void clear() {counter.clear();}
  } membershipCounter;

  // These 3 lists replace the functionality of the "mainQueue" in the
  // SignalQueue ancestor class:
  FixedSizePriorityList<SignalPtr> *lists[3];
  IntronQueueIterator iter, iterEnd;
public:
  IntronQueue(ContentType,
	      int capacity,
	      BOOM::Array1D<SinglePhaseComparator*> &);
  virtual ~IntronQueue();
  void switchComparator(BOOM::Array1D<SinglePhaseComparator*> &);
  virtual void switchToIsochore(Isochore *);
  virtual void resetQueue(Isochore *);
  void purge(int beforePosition);
  virtual void addSignal(SignalPtr);
  virtual void updateHoldingQueue(int position);
  virtual BOOM::Iterator<SignalPtr> &begin();
  virtual BOOM::Iterator<SignalPtr> &end();
  BOOM::Map<SignalPtr ,int>::iterator beginUniq() 
    {return membershipCounter.begin();}
  BOOM::Map<SignalPtr,int>::iterator endUniq() 
    {return membershipCounter.end();}
  virtual int numElements();
  virtual void drop(BOOM::Iterator<SignalPtr>&); 
            // ^-- removes a signal from the queue
  virtual void flushAccumulator();
};


#endif
