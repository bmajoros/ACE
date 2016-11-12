/****************************************************************
 ContentRegions.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/VectorSorter.H"
#include "ContentRegions.H"
using namespace std;
using namespace BOOM;



ContentRegions::ContentRegions(const GffTranscript &transcript,int seqLength)
{
  init(transcript,seqLength);
}



const Vector<ContentRegion> &ContentRegions::asVector()
{
  return regions;
}



void ContentRegions::init(const GffTranscript &transcript,int seqLength)
{
  const int numExons=transcript.numExons();
  if(numExons<1) throw "Empty transcript in ContentRegions::init()";
  regions.push_back(ContentRegion(INTERGENIC,0,transcript.getBegin()));
  for(int i=0 ; i<numExons ; ++i) {
    const GffExon &exon=transcript.getIthExon(i);
    ContentType type;
    if(numExons==1) type=SINGLE_EXON;
    else if(i==0) type=INITIAL_EXON;
    else if(i==numExons-1) type=FINAL_EXON;
    else type=INTERNAL_EXON;
    regions.push_back(ContentRegion(type,exon.getBegin(),exon.getEnd()));
    if(i+1<numExons) {
      const GffExon &nextExon=transcript.getIthExon(i+1);
      regions.push_back(ContentRegion(INTRON,exon.getEnd(),nextExon.getBegin()));
    }
  }
  const int numUTR=transcript.numUTR();
  for(int i=0 ; i<numUTR ; ++i) {
    const GffExon &exon=transcript.getIthUTR(i);
    const Strand strand=exon.getStrand();
    const ContentType type=exonTypeToContentType(exon.getExonType(),strand);
    regions.push_back(ContentRegion(type,exon.getBegin(),exon.getEnd()));
    if(i+1<numUTR) {
      const GffExon &nextExon=transcript.getIthUTR(i+1);
      const Strand nextStrand=nextExon.getStrand();
      const ContentType nextType=exonTypeToContentType(nextExon.getExonType(),
						       nextStrand);
      if(isUTR5(type) && isUTR5(nextType) || 
	 isUTR3(type) && isUTR3(nextType)) {
	ContentType intronType=isUTR5(type) ? UTR5_INTRON : UTR3_INTRON;
	regions.push_back(ContentRegion(intronType,exon.getEnd(),nextExon.getBegin()));
      }
    }
  }
  regions.push_back(ContentRegion(INTERGENIC,transcript.getEnd(),seqLength));
  sort();
}



void ContentRegions::sort()
{
  ContentRegionComparator cmp;
  VectorSorter<ContentRegion> sorter(regions,cmp);
  sorter.sortAscendInPlace();
}



bool ContentRegions::findJunction(int pos,const ContentRegion *&preceding,
				  const ContentRegion *&following) const
{
  for(Vector<ContentRegion>::const_iterator cur=regions.begin(), 
	end=regions.end() ; cur!=end ; ++cur) {
    const ContentRegion *region=&*cur;
    if(region->getInterval().getEnd()==pos) {
      preceding=region;
      ++cur;
      following=&*cur;
      return true;
    }
  }
  return false;
}



ContentRegion *ContentRegions::regionOverlapping(int pos) const
{
  for(Vector<ContentRegion>::const_iterator cur=regions.begin(), 
	end=regions.end() ; cur!=end ; ++cur)
    if((*cur).getInterval().contains(pos)) return &*cur;
  return NULL;
}



void ContentRegions::printOn(ostream &os) const
{
  for(Vector<ContentRegion>::const_iterator cur=regions.begin(), 
	end=regions.end() ; cur!=end ; ++cur) {
    os<<*cur;
    if(cur+1!=end) os<<",";
  }
}



ostream &operator<<(ostream &os,const ContentRegions &r)
{
  r.printOn(os);
  return os;
}


