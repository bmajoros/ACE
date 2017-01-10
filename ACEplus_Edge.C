/****************************************************************
 ACEplus_Edge.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ACEplus_Edge.H"
using namespace std;
using namespace BOOM;

ACEplus_Edge::ACEplus_Edge(const String &substrate,ContentType type,
		     LightVertex *from,LightVertex *to,int begin,int end,
		     Strand s,int ID)
  : LightEdge(substrate,type,from,to,begin,end,s,ID)
{
}



StructureChange &ACEplus_Edge::getChange()
{
  return change;
}



