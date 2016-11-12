/****************************************************************
 SignalQueue.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "SignalQueue.H"
#include <iostream>
#include "genezilla.H"
#include "SignalTypeProperties.H"
#include "Accumulator.H"
#include "BOOM/Exceptions.H"


IgnorantComparator SignalQueue::ignorantComparator;


void SignalQueue::addSignal(SignalPtr signal)
{
  Propagator &propagator=signal->getPropagator(contentType);
  int propagatorPosition=propagator.getPosition();
  if(propagatorPosition <= accumulator.getPosition()) {
    if(propagatorPosition<accumulator.getPosition()) INTERNAL_ERROR;
    flushAccumulator();
    mainQueue.push_back(signal);
  }
  else
    holdingQueue.push_back(signal);
}



void SignalQueue::updateHoldingQueue(int position)
{
  // Checks whether any signals in holdingQueue have been passed by
  // the accumulator.  If so, calls flushAccumulator() and then moves
  // those signals into the main queue.

  int accumulatorPosition=accumulator.getPosition();
  BOOM::List<SignalPtr >::iterator cur=holdingQueue.begin(), 
    end=holdingQueue.end();
  bool first=true;
  for(; cur!=end ;)
    {
      SignalPtr signal=*cur;
      if(signal->getPropagator(contentType).getPosition() <= 
	 accumulatorPosition)
	{
	  if(first) {flushAccumulator();first=false;}
	  BOOM::List<SignalPtr >::iterator graduate=cur;
	  ++cur;
	  holdingQueue.erase(graduate);
	  SignalPtr displaced;
	  if(mainQueue.size()==mainQueue.getCapacity())
	    displaced=*mainQueue.begin();
	  bool inserted=mainQueue.insert(signal);
	  continue;
	}
      ++cur;
    }
}



BOOM::Iterator<SignalPtr> &SignalQueue::begin()
{
  mainIter=mainQueue.begin();
  return mainIter;
}



BOOM::Iterator<SignalPtr> &SignalQueue::end()
{
  mainEnd=mainQueue.end();
  return mainEnd;
}



int SignalQueue::numElements()
{
  return holdingQueue.size() + mainQueue.size();
}



void SignalQueue::drop(BOOM::Iterator<SignalPtr > &victim)
{
  mainQueue.erase(static_cast<FSPLIterator<SignalPtr >&>(victim));
  SignalPtr signal=*victim;
}



void SignalQueue::switchComparator(BOOM::Comparator<SignalPtr> *comp)
{
  mainQueue.changeComparator(comp);
}



void SignalQueue::switchToIsochore(Isochore *)
{
  // do nothing -- overloaded in derived classes
}



void SignalQueue::resetQueue(Isochore *)
{
  mainQueue.clear();
  holdingQueue.clear();

  if(isCoding(contentType))
    accumulator.resetScoresForCoding();       // [0,0,0]
  else
    accumulator.resetScoresForNoncoding();    // [0,-inf,-inf]
}



void SignalQueue::purge(int beforePosition)
{
  for(FixedSizePriorityList<SignalPtr>::iterator cur=mainQueue.begin(),
	next, end=mainQueue.end() ; cur!=end ; ) {
    SignalPtr s=*cur;
    next=cur; ++next;
    if(s->getContextWindowEnd()<=beforePosition) mainQueue.erase(cur);
    cur=next;
  }
}


