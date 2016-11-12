/****************************************************************
 SignalPrinter.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SignalPrinter_H
#define INCL_SignalPrinter_H
#include <iostream>
#include "SignalSensor.H"
using namespace std;
using namespace BOOM;


class SignalPrinter {
public:
  static String print(SignalSensor &sensor,int windowPos,const String &genome);
};


#endif

