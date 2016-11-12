/****************************************************************
 EvidenceFilter.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "EvidenceFilter.H"
using namespace std;
using namespace BOOM;


EvidenceFilter::EvidenceFilter(int m,WigBinary *w,RnaJunctions *j)
  : minSupport(m), wig(w), junctions(j)
{
  // ctor
}



EvidenceFilter::~EvidenceFilter()
{
  delete wig;
  delete junctions;
}



bool EvidenceFilter::intronSupported(int begin,int end)
{
  return junctions->getDepth(begin,end)>=minSupport;
}



bool EvidenceFilter::exonSupported(int begin,int end)
{

  // ### THIS NEEDS TO BE CHANGED TO SOME CONSTANT-TIME ALGORITHM!

  for(int pos=begin ; pos<end ; ++pos)
    if(wig->read(pos)<minSupport) return false;
  return true;
}



bool EvidenceFilter::codingSignalSupported(int begin,int end)
{
  return exonSupported(begin,end);
}



bool EvidenceFilter::spliceOutSupported(int pos)
{
  return junctions->getSpliceOutDepth(pos)>=minSupport;
}



bool EvidenceFilter::spliceInSupported(int pos)
{
  return junctions->getSpliceInDepth(pos)>=minSupport;
}






