/*
 ContentType.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 */
#ifndef INCL_ContentType_H
#define INCL_ContentType_H

#include <iostream>
#include "BOOM/String.H"
#include "BOOM/GffExon.H"
#include "BOOM/Strand.H"
#include "BOOM/Set.H"
#include "SignalType.H"
using namespace BOOM;

/*
  WARNING! Do not change the order of elements in this enum:
  it will break existing model files that have these embedded
  as integers.
 */
enum ContentType
{
  UNKNOWN_CONTENT_FORWARD,
  UNKNOWN_CONTENT_REVERSE,
  INITIAL_EXON,
  INTERNAL_EXON,
  FINAL_EXON,
  SINGLE_EXON,
  INTRON,
  INTERGENIC,
  UTR5_INTRON,
  UTR3_INTRON,
  UTR5_INITIAL,
  UTR5_INTERNAL,
  UTR5_FINAL,
  UTR5_SINGLE, // also serves as general 5' UTR marker
  UTR3_INITIAL,
  UTR3_INTERNAL,
  UTR3_FINAL,
  UTR3_SINGLE, // also serves as general 3' UTR marker
  NEG_INITIAL_EXON,
  NEG_INTERNAL_EXON,
  NEG_FINAL_EXON,
  NEG_SINGLE_EXON,
  NEG_INTRON,
  NEG_UTR5_INTRON,
  NEG_UTR3_INTRON,
  NEG_UTR5_INITIAL,
  NEG_UTR5_INTERNAL,
  NEG_UTR5_FINAL,
  NEG_UTR5_SINGLE,
  NEG_UTR3_INITIAL,
  NEG_UTR3_INTERNAL,
  NEG_UTR3_FINAL,
  NEG_UTR3_SINGLE
};
static const int NumContentTypes=33;

const int FIVE_PRIME_UTR=UTR5_SINGLE;
const int THREE_PRIME_UTR=UTR3_SINGLE;
const int NEG_FIVE_PRIME_UTR=NEG_UTR5_SINGLE;
const int NEG_THREE_PRIME_UTR=NEG_UTR3_SINGLE;

ContentType reverseComplement(ContentType);
ContentType dropStrand(ContentType);
istream &operator>>(istream &,ContentType &);
ostream &operator<<(ostream &,const ContentType &);
Strand getStrand(ContentType);
BOOM::String contentTypeToString(ContentType);
BOOM::String contentTypeNiceString(ContentType);
ContentType stringToContentType(const BOOM::String &);
bool isCoding(ContentType);
bool isIntron(ContentType);
bool isIntergenic(ContentType);
bool isUTR(ContentType);
bool isUTR5(ContentType);
bool isUTR3(ContentType);
bool isUTR5intron(ContentType);
bool isUTR3intron(ContentType);
Set<SignalType> leftSignals(ContentType);
Set<SignalType> rightSignals(ContentType);
ContentType exonTypeToContentType(ExonType,Strand);
ExonType contentTypeToExonType(ContentType);

class ContentTypeInitializer
{
  static ContentTypeInitializer initializer;
  ContentTypeInitializer();
};

#endif
