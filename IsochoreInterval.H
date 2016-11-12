/****************************************************************
 IsochoreInterval.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_IsochoreInterval_H
#define INCL_IsochoreInterval_H
#include "BOOM/Comparator.H"

struct IsochoreInterval
{
  int begin, end;
  int isochoreId; // 0-3
  float GC;
  IsochoreInterval(int b,int e,int id,float gc)
    : begin(b),end(e),isochoreId(id),GC(gc) {}
  IsochoreInterval() : begin(-1) {}
};

class IsochoreIntervalCmp : public BOOM::Comparator<IsochoreInterval>
{
  typedef IsochoreInterval T;
  bool equal(T &a,T &b)   {return a.begin==b.begin;}
  bool greater(T &a,T &b) {return a.begin> b.begin;}
  bool less(T &a,T &b)    {return a.begin< b.begin;}
};

#endif
