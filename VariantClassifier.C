/****************************************************************
 VariantClassifier.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VariantClassifier.H"
using namespace std;
using namespace BOOM;



VariantClassifier::VariantClassifier(const Vector<Variant> &variants,
				     RefAlt refAlt,
				     const GffTranscript &transcript,
				     int chromLen,bool revcomp)
  : CLOSENESS_THRESHOLD(50)
{
  GffTranscript t(transcript);
  if(revcomp) t.reverseComplement(chromLen);
  classify(variants,refAlt,t);
}



const Vector<VariantClassifier::VariantInfo> &
VariantClassifier::getCDSvariants() const
{
  return cdsVariants;
}



const Vector<VariantClassifier::VariantInfo> &
VariantClassifier::getFrameshiftVariants() const
{
  return frameshiftVariants;
}



const Vector<VariantClassifier::VariantInfo> &
VariantClassifier::getSpliceSiteVariants() const
{
  return spliceSiteVariants;
}



const Vector<VariantClassifier::VariantInfo> &
VariantClassifier::getNearSpliceVariants() const
{
  return nearSpliceVariants;
}



const Vector<VariantClassifier::VariantInfo> &VariantClassifier::getAll() const
{
  return all;
}



const Vector<VariantClassifier::VariantInfo> &
VariantClassifier::getUTRvariants() const
{
  return utrVariants;
}



const Vector<VariantClassifier::VariantInfo> &
VariantClassifier::getIntronVariants() const
{
  return intronVariants;
}



void VariantClassifier::classify(const Vector<Variant> &variants,
				 RefAlt refAlt,
				 const GffTranscript &transcript)
{
  int transcriptBegin=transcript.getBegin(), transcriptEnd=transcript.getEnd();
  Vector<Interval> introns;
  transcript.getIntrons(introns);
  for(Vector<Variant>::const_iterator cur=variants.begin(), end=
	variants.end() ; cur!=end ; ++cur) {
    const Variant variant=*cur;
    const int pos=refAlt==REF ? variant.refPos : variant.altPos;
    const int endPos=pos+(refAlt==REF ? variant.alleles[0].length() :
			  variant.alleles[1].length());
    Interval variantInterval(pos,endPos);
    VariantInfo info;
    info.variant=variant;
    info.type=getVariantType(variant);
    if(pos<transcriptBegin || pos>=transcriptEnd)
      { info.elem=INTERGENIC; all.push_back(info); continue; }
    for(Vector<Interval>::const_iterator cur=introns.begin(), end=
	  introns.end() ; cur!=end ; ++cur) {
      const Interval &intron=*cur;
      int iBegin=intron.getBegin(), iEnd=intron.getEnd();
      Interval ss1(iBegin,iBegin+2), ss2(iEnd-2,iEnd);
      if(variantInterval.overlaps(ss1) || variantInterval.overlaps(ss2))
	{ info.elem=SPLICE_SITE; all.push_back(info); continue; }
      if(variantInterval.overlaps(intron))
	{ info.elem=INTRON; all.push_back(info); continue; }
    }
    for(Vector<GffExon*>::iterator cur=transcript.getExons(), end=
	  transcript.getExonsEnd() ; cur!=end ; ++cur) {
      const GffExon *exon=*cur;
      if(exon->overlaps(variantInterval)) {
	info.elem=CDS; all.push_back(info);
	int b=exon->getBegin(), e=exon->getEnd();
	int donorDist=INT_MAX, acceptorDist=INT_MAX;
	if(exon->hasDonor()) donorDist=min(abs(pos-e),abs(endPos-e));
	if(exon->hasAcceptor()) acceptorDist=min(abs(pos-b),endPos-b);
	int dist=min(donorDist,acceptorDist);
	info.distanceToSpliceSite=dist;
	continue;
      }
    }
    for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	  transcript.getUTRend() ; cur!=end ; ++cur) {
      const GffExon *exon=*cur;
      if(exon->overlaps(variantInterval))
	{ info.elem=UTR; all.push_back(info); continue; }
    }
  }

  // Now split them out into individual vectors for easy access
  for(Vector<VariantInfo>::iterator cur=all.begin(), end=all.end() ; 
      cur!=end ; ++cur) {
    const VariantInfo &info=*cur;
    if(info.distanceToSpliceSite>=0 && info.distanceToSpliceSite<=
       CLOSENESS_THRESHOLD)
      nearSpliceVariants.push_back(info);
    switch(info.elem) {
    case UTR: 
      utrVariants.push_back(info); break;
    case CDS:
      cdsVariants.push_back(info);
      if(info.indel() && 
	 abs(int(info.variant.alleles[0].length()-
		 info.variant.alleles[1].length())) % 3 != 0)
	frameshiftVariants.push_back(info);
      break;
    case INTRON:
      intronVariants.push_back(info); // ###
      break;
    case SPLICE_SITE:
      spliceSiteVariants.push_back(info); break;
    case INTERGENIC:
      break;
    }
  }
}



VariantType VariantClassifier::getVariantType(const Variant &variant) const
{
  const int L1=variant.alleles[0].length(), L2=variant.alleles[1].length();
  if(L1==0) return VARIANT_INSERTION;
  if(L2==0) return VARIANT_DELETION;
  if(L1==1 && L2==1) return VARIANT_SNP;
  return VARIANT_COMPLEX;
}



void VariantClassifier::setClosenessThreshold(int t)
{
  CLOSENESS_THRESHOLD=t;
}



Essex::CompositeNode *VariantClassifier::makeVariantsNode() const
{
  if(cdsVariants.size()==0 && spliceSiteVariants.size()==0 &&
     nearSpliceVariants.size()==0 && utrVariants.size()==0 &&
     intronVariants.size()==0) return NULL;
  Essex::CompositeNode *parent=new Essex::CompositeNode("variants");
  addVariants(cdsVariants,"CDS-variants",parent);
  addVariants(frameshiftVariants,"frameshift-variants",parent);
  addVariants(spliceSiteVariants,"splice-site-variants",parent);
  addVariants(nearSpliceVariants,"near-splice-variants",parent);
  addVariants(utrVariants,"UTR-variants",parent);
  addVariants(intronVariants,"intron-variants",parent); // ###
  return parent;
}



void VariantClassifier::addVariants(const Vector<VariantInfo> &variants,
				    const String &tag,
				    Essex::CompositeNode *parent) const
{
  if(variants.empty()) return;
  Essex::CompositeNode *node=new Essex::CompositeNode(tag);
  for(Vector<VariantInfo>::const_iterator cur=variants.begin(), 
	end=variants.end() ; cur!=end ; ++cur) {
    const VariantInfo &info=*cur;
    const Variant &v=info.variant;
    String s=v.id+":"+v.chr+":"+v.refPos+":"+v.altPos+":"+v.alleles[0]
      +":"+v.alleles[1];
    node->append(s);
  }  
  parent->append(node);
}











