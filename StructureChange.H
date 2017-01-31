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

enum StructureChangeType {
  STRUCT_NO_CHANGE,
  STRUCT_EXON_SKIPPING,
  STRUCT_INTRON_RETENTION,
  STRUCT_CRYPTIC_SITE
};

struct StructureChange {
  bool exonSkipping;
  bool intronRetention;
  bool crypticSite;
  bool deNovoSite; // a cryptic site created by an alternate allele
  bool crypticExon; // a new exon created by a de novo splice site
  bool regulatoryChange; // any change brought about by altering splicing
                         // regulatory elements
  StructureChange();
  bool anyChange();
  StructureChange &operator+=(const StructureChange &);
  void printOn(ostream &) const;
};

ostream &operator<<(ostream &,const StructureChange &);

#endif

