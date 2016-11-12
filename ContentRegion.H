/****************************************************************
 ContentRegion.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ContentRegion_H
#define INCL_ContentRegion_H
#include <iostream>
#include "BOOM/Interval.H"
#include "BOOM/Comparator.H"
#include "ContentType.H"
using namespace std;
using namespace BOOM;

class ContentRegion {
public:
  ContentRegion();
  ContentRegion(ContentType,int begin,int end);
  ContentType getType() const;
  const Interval &getInterval() const;
  void printOn(ostream &) const;
protected:
  Interval interval;
  ContentType type;
};

ostream &operator<<(ostream &,const ContentRegion &);

class ContentRegionComparator : public Comparator<ContentRegion> { // compares begin only
public:
  bool equal(ContentRegion &a,ContentRegion &b);
  bool greater(ContentRegion &a,ContentRegion &b);
  bool less(ContentRegion &a,ContentRegion &b);
};

#endif

