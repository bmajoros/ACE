/****************************************************************
 StructureChange.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_StructureChange_H
#define INCL_StructureChange_H
#include <iostream>
using namespace std;
using namespace BOOM;


struct StructureChange {
  bool exonSkipping;
  bool intronRetention;
  bool crypticSite;
  StructureChange() : exonSkipping(false), intronRetention(false), 
		      crypticSite(false) {}
  bool anyChange() { return exonSkipping || intronRetention || crypticSite; }
};

#endif

