/**************************************************************
 NthOrderStringIterator.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "NthOrderStringIterator.H"
#include <iostream>


NthOrderStringIterator::NthOrderStringIterator(int N,Alphabet &alphabet)
  : alphabet(alphabet), N(N), alphabetSize(alphabet.getNumElements())
{
  // ctor

  reset();
}



BOOM::String NthOrderStringIterator::getNextString()
{
  BOOM::String *s=current.toString(alphabet);
  BOOM::String retval=*s;
  delete s;

  advance();

  return retval;
}



bool NthOrderStringIterator::done()
{
  return !hasMore;
}



void NthOrderStringIterator::reset()
{
  Symbol z=0;
  for(int i=0 ; i<N ; ++i)
    current.append(z);
  hasMore=true;
}



void NthOrderStringIterator::advance()
{
  for(int i=N-1 ; i>=0 ; --i)
    {
      Symbol &s=current[i];
      ++s;
      if(s<alphabetSize) return;
      s=0;
    }
  hasMore=false;
}


