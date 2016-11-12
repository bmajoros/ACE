/****************************************************************
 VariantClassifier.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VariantClassifier_H
#define INCL_VariantClassifier_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/GffTranscript.H"
#include "BOOM/Essex.H"
#include "Variant.H"
using namespace std;
using namespace BOOM;

/*
  This class assumes forward-strand genes.
 */

class VariantClassifier {
public:
  enum RefAlt {REF, ALT};
  enum FunctionalElement {
    UTR,
    CDS,
    INTRON,
    SPLICE_SITE,
    INTERGENIC
  };
  struct VariantInfo {
    Variant variant;
    VariantType type; // SNP/INSERTION/DELETION/COMPLEX
    FunctionalElement elem;
    int distanceToSpliceSite;
    VariantInfo() : distanceToSpliceSite(-1) {}
    bool indel() { return type==VARIANT_INSERTION || type==VARIANT_DELETION; }
  };
  VariantClassifier(const Vector<Variant> &,RefAlt,const GffTranscript &,
		    int chromLen,bool revcomp);
  const Vector<VariantInfo> &getCDSvariants() const;
  const Vector<VariantInfo> &getSpliceSiteVariants() const;
  const Vector<VariantInfo> &getNearSpliceVariants() const;
  const Vector<VariantInfo> &getAll() const;
  const Vector<VariantInfo> &getUTRvariants() const;
  const Vector<VariantInfo> &getIntronVariants() const;
  const Vector<VariantInfo> &getFrameshiftVariants() const;
  Essex::CompositeNode *makeVariantsNode() const;
  void setClosenessThreshold(int); // for near-splice-site variants
private:
  int CLOSENESS_THRESHOLD;
  Vector<VariantInfo> cdsVariants, spliceSiteVariants, nearSpliceVariants;
  Vector<VariantInfo> all, utrVariants, intronVariants, frameshiftVariants;
  void classify(const Vector<Variant> &,RefAlt,const GffTranscript &);
  VariantType getVariantType(const Variant &) const;
  void addVariants(const Vector<VariantInfo> &,const String &tag,
		   Essex::CompositeNode *parent) const;
};

#endif

