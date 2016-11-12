/****************************************************************
 Variant.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Variant.H"
using namespace std;
using namespace BOOM;


ostream &operator<<(ostream &os,const Variant &v) 
{
  v.printOn(os);
  return os;
}




