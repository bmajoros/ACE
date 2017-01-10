/****************************************************************
 ACEplus_Edge.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ACEplus_Edge_H
#define INCL_ACEplus_Edge_H
#include <iostream>
#include "LightEdge.H"
#include "StructureChange.H"
using namespace std;
using namespace BOOM;


class ACEplus_Edge : public LightEdge {
public:
  ACEplus_Edge(const String &substrate,ContentType,LightVertex *,LightVertex *,
	    int begin,int end,Strand,int ID);
  StructureChange &getChange();
protected:
  StructureChange change;
};

#endif

