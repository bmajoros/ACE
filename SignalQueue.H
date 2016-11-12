/****************************************************************
 SignalQueue.H Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalQueue_H
#define INCL_SignalQueue_H
#include "BOOM/FixedSizePriorityList.H"
#include "BOOM/Constants.H"
#include "ContentType.H"
#include "Signal.H"
#include "Accumulator.H"
#include "SignalComparator.H"
#include "Isochore.H"
#include "DiscreteDistribution.H"

/************************************************************************
 Each SignalQueue contains only signals which are followed by the same
 type of content in a valid gene model; thus, each queue has its own
 content sensor for evaluating the content following all the signals
 in the queue [*** that's not quite true anymore... ***].

 The accumulator moves along the sequence accumulating scores from
 individual bases as they are evaluated by this queue's content sensor.
 At certain key times the contents of the accumulator are added to the
 propagators of all the signals in this queue (and the accumulator score is
 then reset to zero).  Thus we maintain the invariant that a signal's
 propagator, appropriately added to the accumulator (i.e., accounting
 for phase/frame), gives the inductive score of the best gene model
 up to this point, passing through that signal.
 
 The main queue contains all signals of a given type that have been
 completely passed by this queue's accumulator; they therefore are
 eligible to be linked to any subsequent signal that is found.
 
 The holdingQueue is for those signals which have been very recently
 added to the queue, so that their propagators (by virtue of having
 been initialized to include the scores in the entire context window)
 are currently ahead of this queue's accumulator.  When that is no
 longer true they are moved to the main queue and are then eligible
 to be the predecessor of a subsequent signal.
***********************************************************************/

class SignalQueue
{
  static IgnorantComparator ignorantComparator;
protected:
  ContentType contentType;
  Accumulator accumulator;
  FixedSizePriorityList<SignalPtr> mainQueue;
  BOOM::List<SignalPtr> holdingQueue;
  FSPLIterator<SignalPtr> mainIter, mainEnd;
  BOOM::List<SignalPtr>::iterator holdingIter,holdingEnd;
public:
  inline SignalQueue(ContentType,
		     int capacity=LARGEST_INTEGER,
		     BOOM::Comparator<SignalPtr> *comp=&ignorantComparator);
  virtual ~SignalQueue() 
    {/* signals are deleted by GarbageCollector, not here! */}
  FixedSizePriorityList<SignalPtr> &getMainQueue() {return mainQueue;}
  void switchComparator(BOOM::Comparator<SignalPtr> *);
  virtual void switchToIsochore(Isochore *);
  virtual void resetQueue(Isochore *);
  virtual void purge(int beforePosition);
  inline bool precedesExon();
  inline bool precedesIntron();
  inline bool precedesIntergenic();
  inline bool precedesUTR();
  inline Strand getStrand();
  inline ContentType getContentType();
  virtual void addSignal(SignalPtr );
  inline Accumulator &getAccumulator();
  virtual BOOM::Iterator<SignalPtr> &begin();
  virtual BOOM::Iterator<SignalPtr> &end();
  inline BOOM::List<SignalPtr>::iterator &getHoldingQueue();
  inline BOOM::List<SignalPtr>::iterator &holdingQueueEnd();
  virtual int numElements();
  inline bool isEmpty();
  virtual void drop(BOOM::Iterator<SignalPtr> &); // removes signal from queue
  inline void dropFromHoldingQueue(BOOM::List<SignalPtr>::iterator &);

  inline void addToAccumulator(double scorePhase0,double scorePhase1,
			double scorePhase2,int positionOnSubstrate);
  inline void addToAccumulator(double score);
  // ^--Adds values to accumulator according to absolute (substrate) frame,
  //    advances accumulator's position, and then calls updateHoldingQueue()

  virtual void updateHoldingQueue(int position);
  // ^--Checks whether any signals in holdingQueue have been passed by
  //    the accumulator.  If so, calls flushAccumulator() and then moves
  //    those signals into the main queue.

  virtual void flushAccumulator();
  // ^--Adds accumulator contents to propagators of all signals in the
  //    main queue (after adjusting for frame/phase) and then resets
  //    accumulator scores to 0.

  void resetPosition() {accumulator.resetPosition();}
};



inline SignalQueue::SignalQueue(ContentType contentType,
				int capacity,
				BOOM::Comparator<SignalPtr> *comp)
  : accumulator(-1),
    mainQueue(capacity,comp),
    contentType(contentType)
{
}



inline bool SignalQueue::precedesExon()
{
  return isCoding(contentType);
}



inline bool SignalQueue::precedesIntron()
{
  return isIntron(contentType);
}



inline bool SignalQueue::precedesIntergenic()
{
  return isIntergenic(contentType);
}



inline bool SignalQueue::precedesUTR()
{
  return isUTR(contentType);
}



inline ContentType SignalQueue::getContentType()
{
  return contentType;
}



inline Accumulator &SignalQueue::getAccumulator()
{
  return accumulator;
}



inline BOOM::List<SignalPtr >::iterator &SignalQueue::getHoldingQueue()
{
  holdingIter=holdingQueue.begin();
  return holdingIter;
}



inline BOOM::List<SignalPtr >::iterator &SignalQueue::holdingQueueEnd()
{
  holdingEnd=holdingQueue.end();
  return holdingEnd;
}



inline void SignalQueue::addToAccumulator(double scorePhase0,
					  double scorePhase1,
					  double scorePhase2,int position)
{
  // CODING VERSION

  // Adds values to accumulator according to absolute (substrate) frame,
  // advances accumulator's position, and then calls updateHoldingQueue().

  int absoluteFrame=position%3;//### I think the % is unnecessary...
  switch(getStrand())
    {
    case FORWARD_STRAND:
      accumulator[posmod(0-absoluteFrame)]+=scorePhase0;
      accumulator[posmod(1-absoluteFrame)]+=scorePhase1;
      accumulator[posmod(2-absoluteFrame)]+=scorePhase2;      
      break;

    case REVERSE_STRAND:
      accumulator[posmod(0+absoluteFrame)]+=scorePhase0;
      accumulator[posmod(1+absoluteFrame)]+=scorePhase1;
      accumulator[posmod(2+absoluteFrame)]+=scorePhase2;      
      break;
    }
  accumulator.incrementPosition();
  updateHoldingQueue(accumulator.getPosition());
}



inline Strand SignalQueue::getStrand()
{
  return ::getStrand(contentType);
}



inline void SignalQueue::addToAccumulator(double score)
{
  // NONCODING VERSION

  // Adds values to accumulator, advances accumulator's position, and 
  // then calls updateHoldingQueue().

  accumulator.addFrameZeroAndAdvance(score);
  updateHoldingQueue(accumulator.getPosition());
}



inline void SignalQueue::flushAccumulator()
{
  // Adds accumulator contents to propagators of all signals in the
  // main queue (after adjusting for frame/phase) and then resets
  // accumulator scores to 0.

  // Update the signal propagators:
  Strand strand=getStrand();
  bool isCoding=precedesExon();
  FixedSizePriorityList<SignalPtr >::iterator cur=mainQueue.begin(), 
    end=mainQueue.end();
  for(; cur!=end ; ++cur)
    {
      SignalPtr signal=*cur;
      Propagator &propagator=signal->getPropagator(contentType);
      if(signal->getSignalType()==ATG) cout<<"XXX "<<accumulator<<" "<<propagator<<" => ";
      propagator.update(accumulator,strand,isCoding,signal);      
      if(signal->getSignalType()==ATG) cout<<propagator<<endl;
    }

  // Reset accumulator scores:
  if(isCoding)
    accumulator.resetScoresForCoding();       // [0,0,0]
  else
    accumulator.resetScoresForNoncoding();    // [0,-inf,-inf]
}



inline void SignalQueue::dropFromHoldingQueue
  (BOOM::List<SignalPtr>::iterator &victim)
{
  holdingQueue.erase(victim);
  SignalPtr signal=*victim;
}



inline bool SignalQueue::isEmpty()
{
  return numElements()==0;
}



#endif
