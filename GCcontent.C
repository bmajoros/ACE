/****************************************************************
 GCcontent.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "GCcontent.H"
using namespace std;
using namespace BOOM;


float GCcontent::get(const String &s)
{
  const int A=s.count('A')+s.count('a');
  const int C=s.count('C')+s.count('c');
  const int G=s.count('G')+s.count('g');
  const int T=s.count('T')+s.count('t');
  const int total=A+C+G+T;
  const float gc=(G+C)/float(total);
  return gc;
}


