/****************************************************************
 EvidenceFilter.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_EvidenceFilter_H
#define INCL_EvidenceFilter_H
#include <iostream>
#include "BOOM/WigBinary.H"
#include "RnaJunctions.H"
using namespace std;
using namespace BOOM;

class EvidenceFilter {
public:
  EvidenceFilter(int minSupport,WigBinary*,RnaJunctions*);
  virtual ~EvidenceFilter();
  virtual bool intronSupported(int begin,int end); // 0-based, half open
  virtual bool exonSupported(int begin,int end);   // 0-based, half open
  virtual bool codingSignalSupported(int begin,int end); // ditto
  virtual bool spliceOutSupported(int pos);
  virtual bool spliceInSupported(int pos);
protected:
  int minSupport;
  WigBinary *wig;
  RnaJunctions *junctions;
};

#endif

