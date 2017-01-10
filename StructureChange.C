/****************************************************************
 StructureChange.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/Exceptions.H"
#include "StructureChange.H"
using namespace std;
using namespace BOOM;


StructureChange::StructureChange()
  : exonSkipping(false), intronRetention(false), crypticSite(false)
{
  // ctor
}



bool StructureChange::anyChange()
{
  return exonSkipping || intronRetention || crypticSite;
}



StructureChange &StructureChange::operator+=(const StructureChange &other)
{
  crypticSite=crypticSite || other.crypticSite;
  exonSkipping=exonSkipping || other.exonSkipping;
  intronRetention=intronRetention || other.intronRetention;
}



void StructureChange::printOn(ostream &os) const
{
  int count=0;
  if(crypticSite) { os<<"cryptic-splicing"; ++count; }
  if(exonSkipping) {
    if(count>0) os<<",";
    os<<"exon-skipping";
    ++count; }
  if(intronRetention) {
    if(count>0) os<<",";
    os<<"intron-retention";
    ++count; }
}



ostream &operator<<(ostream &os,const StructureChange &c)
{
  c.printOn(os);
  return os;
}



