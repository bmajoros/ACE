/****************************************************************
 RnaJunctions.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "RnaJunctions.H"
#include "BOOM/VectorSorter.H"
using namespace std;
using namespace BOOM;

RnaJunctions::RnaJunctions(const String &filename)
{
  if(!filename.isEmpty()) load(filename);
}



bool RnaJunctions::load(const String &filename)
{
  File f(filename);
  long fileSize=f.getSize();
  int recSize=3*sizeof(int)+sizeof(char);
  int numRecords=fileSize/recSize;
  junctions.resize(numRecords);
  for(int i=0 ; i<numRecords ; ++i) 
    if(!junctions[i].read(f)) return false;
  sortEnds();
  return true;
}



void RnaJunctions::sortEnds()
{
  JunctionEndComparator cmp;
  VectorSorter<RnaJunction> sorter(junctions,cmp);
  sorter.sortAscendByIndex(sortedByEnd);
}



int RnaJunctions::getNumJunctions() const
{
  return junctions.size();
}



const RnaJunction &RnaJunctions::operator[](int i) const
{
  return junctions[i];
}



float RnaJunctions::getDepth(int begin,int end) const
{
  // Start with a binary-search like procedure to find the last
  // record having the right "begin" value
  const int N=junctions.size();
  if(N==0) return 0;
  int b=0, e=N;
  while(b+1<e) {
    int m=(b+e)/2;
    const RnaJunction &J=junctions[m];
    if(begin<J.getBegin()) e=m;
    else b=m;
  }
  // Now if there are any records with the appropriate "begin" value,
  // they must lie immediately to the left of e
  for(int i=e-1 ; i>=0  ; --i) {
    const RnaJunction &J=junctions[i];
    if(J.getBegin()!=begin) break;
    if(J.getEnd()==end) return J.getDepth();
  }
  return 0;
}





float RnaJunctions::getSpliceOutDepth(int pos) const // leftmost intron pos
{
  // Start with a binary-search like procedure to find the last
  // record having the right "begin" value
  float sum=0;
  const int N=junctions.size();
  if(N==0) return 0;
  int b=0, e=N;
  while(b+1<e) {
    int m=(b+e)/2;
    const RnaJunction &J=junctions[m];
    if(pos<J.getBegin()) e=m;
    else b=m;
  }
  // Now if there are any records with the appropriate "begin" value,
  // they must lie immediately to the left of e
  for(int i=e-1 ; i>=0  ; --i) {
    const RnaJunction &J=junctions[i];
    if(J.getBegin()!=pos) break;
    sum+=J.getDepth();
  }
  return sum;
}



float RnaJunctions::getSpliceInDepth(int pos) const // rightmost intron pos + 1
{
  // Start with a binary-search like procedure to find the last
  // record having the right "end" value
  float sum=0;
  const int N=junctions.size();
  if(N==0) return 0;
  int b=0, e=N;
  while(b+1<e) {
    int m=(b+e)/2;
    const RnaJunction &J=junctions[sortedByEnd[m]];
    if(pos<J.getEnd()) e=m;
    else b=m;
  }
  // Now if there are any records with the appropriate "end" value,
  // they must lie immediately to the left of e
  for(int i=e-1 ; i>=0  ; --i) {
    const RnaJunction &J=junctions[sortedByEnd[i]];
    if(J.getEnd()!=pos) break;
    sum+=J.getDepth();
  }
  return sum;
}



float RnaJunctions::getMaxDepth() const
{
  float max=0;
  Vector<RnaJunction>::const_iterator cur=junctions.begin(), 
    end=junctions.end();
  for(; cur!=end ; ++cur) {
    float depth=(*cur).getDepth();
    if(depth>max) max=depth;
  }
  return max;
}


void RnaJunctions::logify()
{
  Vector<RnaJunction>::iterator cur=junctions.begin(), end=junctions.end();
  for(; cur!=end ; ++cur) (*cur).logify();
}

