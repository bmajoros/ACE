/****************************************************************
 NMD.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_NMD_H
#define INCL_NMD_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/GffTranscript.H"
using namespace std;
using namespace BOOM;

enum NMD_TYPE {
  NMD_NONE,       // nothing wrong
  NMD_NMD,        // premature stop codon & nonsense-medicated decay
  NMD_TRUNCATION, // premature stop codon, truncated protein
  NMD_NO_STOP,    // no stop codon
  NMD_NO_START    // no start codon
};
typedef NMD_TYPE ProteinFate;
ostream &operator<<(ostream &,NMD_TYPE);

class NMD {
public:
  NMD(int dist=50);
  NMD_TYPE predict(GffTranscript &,const String &genome,int &distance);
  void setDistParm(int d) { distParm=d; }
protected:
  int distParm;
};

#endif

