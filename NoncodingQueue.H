/****************************************************************
 NoncodingQueue.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_NoncodingQueue_H
#define INCL_NoncodingQueue_H

#include "IntronQueue.H"

/*
  This is a special type of SignalQueue in which there can be multiple
  (overlapping) signals in the holding queue, but there can be only
  one signal in the main queue.  When a signal is ready to graduate
  from the holding queue to the main queue we see if it is better than
  the one in the main queue (i.e., has a better inductive score), and 
  if so, move it there; otherwise we discard it.
 */

class NoncodingQueue : public SignalQueue
{
public:
  NoncodingQueue(ContentType,
		 int capacity,
		 NoncodingComparator *);
  virtual void updateHoldingQueue(int position);
  virtual void switchToIsochore(Isochore *);
  virtual void resetQueue(Isochore *);
};


#endif
