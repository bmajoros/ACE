/****************************************************************
 RnaJunction.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "RnaJunction.H"
using namespace std;
using namespace BOOM;

RnaJunction::RnaJunction()
{
  // ctor
}



bool RnaJunction::read(File &f)
{
  begin=f.readInt();
  if(f.eof()) return false;
  end=f.readInt();
  depth=f.readFloat();
  strand=f.readChar();
  return true;
}



bool RnaJunction::write(File &f)
{
  f.write(begin);
  f.write(end);
  f.write(depth);
  f.write(strand);
  return true;
}



int RnaJunction::getBegin()
{
  return begin;
}



int RnaJunction::getEnd()
{
  return end;
}



float RnaJunction::getDepth()
{
  return depth;
}



Strand RnaJunction::getStrand()
{
  return strand;
}



void RnaJunction::printOn(ostream &os) const
{
  //os<<begin<<"-"<<end<<"("<<depth<<":"<<strand<<")";
  os<<begin<<"\t"<<end<<"\t"<<depth<<"\t"<<strand;
}



ostream &operator<<(ostream &os,const RnaJunction &J)
{
  J.printOn(os);
  return os;
}



