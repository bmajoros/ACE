/****************************************************************
 SignalStream.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SignalStream.H"
#include "BOOM/VectorSorter.H"
using namespace std;
using namespace BOOM;



SignalStream::SignalStream()
  : currentIndex(0)
{
  // ctor
}



SignalStream::~SignalStream()
{
  // don't delete the signals: the garbage collector does that
}



void SignalStream::add(Signal *s)
{
  signals.push_back(s);
}



void SignalStream::sort()
{
  SignalPosComparator cmp;
  VectorSorter<Signal*> sorter(signals,cmp);
  sorter.sortAscendInPlace();
}



Signal *SignalStream::detect(int position)
{
  int N=signals.size();
  if(currentIndex>=N) return NULL;
  const sigPos=signals[currentIndex]->getContextWindowPosition();
  if(position==sigPos) return signals[currentIndex++];
  if(position>sigPos) {
    cout<<position<<" > "<<sigPos<<" next sig="<<*signals[currentIndex]<<endl;
    INTERNAL_ERROR; // ### DEBUGGING
  }
  return NULL;
}



void SignalStream::reset()
{
  currentIndex=0;
}



void SignalStream::deduplicate()
{
  int N=signals.size();
  for(int i=0 ; i+1<N ; ) {
    Signal *thisOne=signals[i], *nextOne=signals[i+1];
    if(thisOne->getContextWindowPosition()==
       nextOne->getContextWindowPosition() &&
       thisOne->getSignalType()==
       nextOne->getSignalType() &&
       thisOne->getStrand()==
       nextOne->getStrand())
      { signals.cut(i); --N; }
    else ++i;
  }
}



void SignalStream::printOn(ostream &os) const
{
  //  Vector<Signal*> signals;
  for(Vector<Signal*>::const_iterator cur=signals.begin(), 
	end=signals.end() ; cur!=end; ++cur)
    os<<**cur<<" ";
}



ostream &operator<<(ostream &os,const SignalStream &s)
{
  s.printOn(os);
  return os;
}




