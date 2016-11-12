/****************************************************************
 IntronQueue.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "BOOM/Constants.H"
#include "IntronQueue.H"
#include <iostream>
#include "genezilla.H"


IntronQueue::IntronQueue(ContentType contentType,
			 int capacity,
			 BOOM::Array1D<SinglePhaseComparator*> &comp)
  : SignalQueue(contentType)
{
  for(int i=0 ; i<3 ; ++i)
    lists[i]=new FixedSizePriorityList<SignalPtr>(capacity,comp[i]);
  accumulator[1]=accumulator[2]=NEGATIVE_INFINITY;
}



IntronQueue::~IntronQueue()
{
  delete lists[0];
  delete lists[1];
  delete lists[2];

  // do not delete signals here!  they are deleted elsewhere!
}



void IntronQueue::switchComparator(BOOM::Array1D<SinglePhaseComparator*> &cmp)
{
  for(int i=0 ; i<3 ; ++i)
    lists[i]->changeComparator(cmp[i]);
}



void IntronQueue::switchToIsochore(Isochore *isochore)
{
  ContentType contentType=getContentType();
  BOOM::Array1D<SinglePhaseComparator*> *cmp=
    isochore->intronComparators[contentType];
  switchComparator(*cmp);
}



void IntronQueue::resetQueue(Isochore *isochore)
{
  mainQueue.clear();
  holdingQueue.clear();
  for(int i=0 ; i<3 ; ++i)
    lists[i]->clear();
  membershipCounter.clear();

  ContentType contentType=getContentType();
  BOOM::Array1D<SinglePhaseComparator*> *cmp=
    isochore->intronComparators[contentType];
  switchComparator(*cmp);

  accumulator.resetScoresForNoncoding();       // [0,-inf,-inf]
}



void IntronQueue::purge(int beforePosition)
{
  for(int i=0 ; i<3 ; ++i) {
    FixedSizePriorityList<SignalPtr> &q=*lists[i];
    for(FixedSizePriorityList<SignalPtr>::iterator cur=q.begin(), next,
	  end=q.end() ; cur!=end ; ) {
      SignalPtr s=*cur;
      next=cur; ++next;
      if(s->getContextWindowEnd()<=beforePosition) {
	membershipCounter.decrement(s);
	mainQueue.erase(cur);
      }
      cur=next;
    }
  }
}



void IntronQueue::addSignal(SignalPtr signal)
{
  Propagator &propagator=signal->getPropagator(contentType);
  int propagatorPosition=propagator.getPosition();
  if(propagatorPosition <= accumulator.getPosition())
                     // ^-- only true for the "left terminus" signals
    for(int i=0 ; i<3 ; ++i)
      {
	lists[i]->push_back(signal);
	membershipCounter.increment(signal);
      }
  else
    holdingQueue.push_back(signal);
}



void IntronQueue::updateHoldingQueue(int position)
{
  // Checks whether any signals in holdingQueue have been passed by
  // the accumulator.  If so, calls flushAccumulator() and then moves
  // those signals into the main queue.

  int accumulatorPosition=accumulator.getPosition();

  // Iterate through the holding queue
  bool changes=false;
  BOOM::List<SignalPtr >::iterator cur=holdingQueue.begin(), 
    end=holdingQueue.end();
  bool first=true;
  while(cur!=end)
    {
      SignalPtr signal=*cur;
      Propagator &newProp=signal->getPropagator(contentType);
      if(newProp.getPosition() <= accumulatorPosition)
	{
	  changes=true;
	  // This signal is ready to be shifted into the main queue.

	  // Make sure any signals already in the main queue are propagated
	  // up to this point before proceding with the promotion (so that
	  // inductive scores will be comparable):
	  if(first) {flushAccumulator();first=false;}

	  // Remove it from the holding queue
	  BOOM::List<SignalPtr >::iterator graduate=cur;
	  ++cur;
	  holdingQueue.erase(graduate);

	  // Add it to the three phase-specific priority queues
	  for(int phase=0 ; phase<3 ; ++phase)
	    {
	      double newScore=newProp[phase];
	      if(isinf(newScore)) continue;
	      FixedSizePriorityList<SignalPtr > *q=lists[phase];
	      SignalPtr displaced;
	      if(q->size()==q->getCapacity()) displaced=*q->begin();
	      bool inserted=q->insert(signal);
                                //may reject it if score is too low
	      if(inserted) 
		{
		  membershipCounter.increment(signal);
		  if(displaced) 
		      membershipCounter.decrement(displaced);
		}
	    }
	}
      else ++cur;
    }
}



BOOM::Iterator<SignalPtr > &IntronQueue::begin()
{
  iter=IntronQueueIterator(lists,false); // false = not past-the-end
  return iter;
}



BOOM::Iterator<SignalPtr > &IntronQueue::end()
{
  iterEnd=IntronQueueIterator(lists,true); // true = past-the-end
  return iterEnd;
}



int IntronQueue::numElements()
{
  int n=0;
  for(int i=0 ; i<3 ; ++i)
    n+=lists[i]->size();
  return n;
}



void IntronQueue::drop(BOOM::Iterator<SignalPtr > &victim)
{
  SignalPtr signal=*victim;
  membershipCounter.decrement(signal);
  IntronQueueIterator &iqi=static_cast<IntronQueueIterator&>(victim);
  lists[iqi.getIndex()]->erase(iqi.getNativeIterator());
}



void IntronQueue::flushAccumulator()
{
  // Adds accumulator contents to propagators of all signals in the
  // main queue and then resets accumulator scores to 0.

  // Update the signal propagators:
  Strand strand=getStrand();

  BOOM::Map<SignalPtr ,int>::iterator cur=membershipCounter.begin(),
    end=membershipCounter.end();
  for(; cur!=end ; ++cur)
    {
      SignalPtr signal=(*cur).first;
      Propagator &propagator=signal->getPropagator(contentType);
      propagator.update(accumulator,strand,false,signal);      
    }

  // Reset accumulator scores:
  accumulator.resetScoresForNoncoding();       // [0,-inf,-inf]
}



//==============================================================
//               IntronQueueIterator methods
//==============================================================

IntronQueueIterator::IntronQueueIterator(FixedSizePriorityList<SignalPtr > 
					 *lists[3],
					 bool pastTheEnd)
{
  if(pastTheEnd) 
    {
      index=2;
      iterators[2]=lists[2]->end();
    }
  else
    {
      index=0;
      for(int i=0 ; i<3 ; ++i)
	{
	  iterators[i]=lists[i]->begin();
	  ends[i]=lists[i]->end();
	}
    }
}



void IntronQueueIterator::operator++()
{
  increment();
}



void IntronQueueIterator::increment()  
{
  FixedSizePriorityList<SignalPtr>::iterator &iter=iterators[index];
  ++iter;
  if(iter==ends[index] && index<2) ++index;
}



void IntronQueueIterator::operator++(int)
{
  increment();
}



SignalPtr &IntronQueueIterator::operator*()
{
  return *iterators[index];
}



bool IntronQueueIterator::operator==(const BOOM::Iterator<SignalPtr> &other) const
{
  const IntronQueueIterator &iqiOther=
    dynamic_cast<const IntronQueueIterator&>(other);
  return iterators[index]==iqiOther.iterators[index];
}



FixedSizePriorityList<SignalPtr>::iterator 
   IntronQueueIterator::getNativeIterator()
  const 
{
  return iterators[index];
}



int IntronQueueIterator::getIndex() const
{
  return index;
}



BOOM::Iterator<SignalPtr> &IntronQueueIterator::clone()
{
  return *new IntronQueueIterator(*this);
}




