/****************************************************************
 NoncodingQueue.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "NoncodingQueue.H"
#include <iostream>
#include "BOOM/Constants.H"
#include "genezilla.H"


NoncodingQueue::NoncodingQueue(ContentType contentType,
			       int capacity,
			       NoncodingComparator *comp)
  : SignalQueue(contentType,capacity,comp)
{
  // base-class ctor initializes all 3 entries to zero
  accumulator[1]=accumulator[2]=NEGATIVE_INFINITY;
}



void NoncodingQueue::updateHoldingQueue(int position)
{
  // Checks whether any signals in holdingQueue have been passed by
  // the accumulator.  If so, calls flushAccumulator() and then moves
  // those signals into the main queue.

  int accumulatorPosition=accumulator.getPosition();

  // Iterate through the holding queue
  BOOM::List<SignalPtr >::iterator cur=holdingQueue.begin(), 
    end=holdingQueue.end();
  bool first=true;
  while(cur!=end)
    {
      SignalPtr signal=*cur;
      Propagator &newProp=signal->getPropagator(contentType);
      if(newProp.getPosition() <= accumulatorPosition)
	{
	  // This signal is ready to be shifted into the main queue.

	  // Make sure any signal already in the main queue is propagated
	  // up to this point before proceding with the promotion (so that
	  // inductive scores will be comparable):
	  if(first) {flushAccumulator();first=false;}

	  // Remove it from the holding queue
	  BOOM::List<SignalPtr >::iterator graduate=cur;
	  ++cur;
	  holdingQueue.erase(graduate);
	  /*
	  SignalPtr displaced;
	  if(mainQueue.size()==mainQueue.getCapacity()) 
	    displaced=*mainQueue.begin();
	  */

	  // Add it to the main queue
	  bool inserted=mainQueue.insert(signal);
	}
      else ++cur;
    }
}



void NoncodingQueue::switchToIsochore(Isochore *isochore)
{
  switchComparator(isochore->noncodingComparators[getContentType()]);
}



void NoncodingQueue::resetQueue(Isochore *isochore)
{
  mainQueue.clear();
  holdingQueue.clear();
  switchComparator(isochore->noncodingComparators[getContentType()]);

  accumulator.resetScoresForNoncoding();       // [0,-inf,-inf]
}


