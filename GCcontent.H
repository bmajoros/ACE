/****************************************************************
 GCcontent.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_GCcontent_H
#define INCL_GCcontent_H
#include <iostream>
#include "BOOM/String.H"
using namespace std;
using namespace BOOM;

class GCcontent {
public:
  static float get(const String &); // returns value between 0 and 1
};

#endif

