/*
 ContentType.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 */
#include "ContentType.H"
#include "BOOM/Map.H"


static BOOM::Map<BOOM::String,ContentType> contentTypeMap;
static BOOM::Map<BOOM::String,ContentType> longContentTypeMap;
ContentTypeInitializer ContentTypeInitializer::initializer;


ContentTypeInitializer::ContentTypeInitializer()
{
  contentTypeMap["INTERGENIC"]=      INTERGENIC;
  contentTypeMap["INITIAL-EXON"]=    INITIAL_EXON;
  contentTypeMap["INTERNAL-EXON"]=   INTERNAL_EXON;
  contentTypeMap["FINAL-EXON"]=      FINAL_EXON;
  contentTypeMap["SINGLE-EXON"]=     SINGLE_EXON;
  contentTypeMap["INTRON"]=          INTRON;
  contentTypeMap["UTR5-INTRON"]=     UTR5_INTRON;
  contentTypeMap["UTR3-INTRON"]=     UTR3_INTRON;
  contentTypeMap["INITIAL-UTR5"]=    UTR5_INITIAL;
  contentTypeMap["INTERNAL-UTR5"]=   UTR5_INTERNAL;
  contentTypeMap["FINAL-UTR5"]=      UTR5_FINAL;
  contentTypeMap["SINGLE-UTR5"]=     UTR5_SINGLE;
  contentTypeMap["INITIAL-UTR3"]=    UTR3_INITIAL;
  contentTypeMap["INTERNAL-UTR3"]=   UTR3_INTERNAL;
  contentTypeMap["FINAL-UTR3"]=      UTR3_FINAL;
  contentTypeMap["SINGLE-UTR3"]=     UTR3_SINGLE;

  contentTypeMap["UTR5-INITIAL"]=    UTR5_INITIAL;
  contentTypeMap["UTR5-INTERNAL"]=   UTR5_INTERNAL;
  contentTypeMap["UTR5-FINAL"]=      UTR5_FINAL;
  contentTypeMap["UTR5-SINGLE"]=     UTR5_SINGLE;
  contentTypeMap["UTR3-INITIAL"]=    UTR3_INITIAL;
  contentTypeMap["UTR3-INTERNAL"]=   UTR3_INTERNAL;
  contentTypeMap["UTR3-FINAL"]=      UTR3_FINAL;
  contentTypeMap["UTR3-SINGLE"]=     UTR3_SINGLE;

  contentTypeMap["NEG-INITIAL-UTR5"]=    NEG_UTR5_INITIAL;
  contentTypeMap["NEG-INTERNAL-UTR5"]=   NEG_UTR5_INTERNAL;
  contentTypeMap["NEG-FINAL-UTR5"]=      NEG_UTR5_FINAL;
  contentTypeMap["NEG-SINGLE-UTR5"]=     NEG_UTR5_SINGLE;
  contentTypeMap["NEG-INITIAL-UTR3"]=    NEG_UTR3_INITIAL;
  contentTypeMap["NEG-INTERNAL-UTR3"]=   NEG_UTR3_INTERNAL;
  contentTypeMap["NEG-FINAL-UTR3"]=      NEG_UTR3_FINAL;
  contentTypeMap["NEG-SINGLE-UTR3"]=     NEG_UTR3_SINGLE;

  contentTypeMap["NEG-UTR5-INITIAL"]=    NEG_UTR5_INITIAL;
  contentTypeMap["NEG-UTR5-INTERNAL"]=   NEG_UTR5_INTERNAL;
  contentTypeMap["NEG-UTR5-FINAL"]=      NEG_UTR5_FINAL;
  contentTypeMap["NEG-UTR5-SINGLE"]=     NEG_UTR5_SINGLE;
  contentTypeMap["NEG-UTR3-INITIAL"]=    NEG_UTR3_INITIAL;
  contentTypeMap["NEG-UTR3-INTERNAL"]=   NEG_UTR3_INTERNAL;
  contentTypeMap["NEG-UTR3-FINAL"]=      NEG_UTR3_FINAL;
  contentTypeMap["NEG-UTR3-SINGLE"]=     NEG_UTR3_SINGLE;

  contentTypeMap["NEG-INITIAL-EXON"]=    NEG_INITIAL_EXON;
  contentTypeMap["NEG-INTERNAL-EXON"]=   NEG_INTERNAL_EXON;
  contentTypeMap["NEG-FINAL-EXON"]=      NEG_FINAL_EXON;
  contentTypeMap["NEG-SINGLE-EXON"]=     NEG_SINGLE_EXON;
  contentTypeMap["NEG-INTRON"]=          NEG_INTRON;
  contentTypeMap["NEG-UTR5-INTRON"]=     NEG_UTR5_INTRON;
  contentTypeMap["NEG-UTR3-INTRON"]=     NEG_UTR3_INTRON;
  contentTypeMap["UNKNOWN-CONTENT-FORWARD"]=UNKNOWN_CONTENT_FORWARD;
  contentTypeMap["UNKNOWN-CONTENT-REVERSE"]=UNKNOWN_CONTENT_REVERSE;

  longContentTypeMap["initial-exon"]=    INITIAL_EXON;
  longContentTypeMap["internal-exon"]=   INTERNAL_EXON;
  longContentTypeMap["final-exon"]=      FINAL_EXON;
  longContentTypeMap["single-exon"]=     SINGLE_EXON;
  longContentTypeMap["intron"]=          INTRON;
  longContentTypeMap["UTR5-intron"]=     UTR5_INTRON;
  longContentTypeMap["UTR3-intron"]=     UTR3_INTRON;
  longContentTypeMap["intergenic"]=      INTERGENIC;
  longContentTypeMap["initial-UTR5"]=    UTR5_INITIAL;
  longContentTypeMap["internal-UTR5"]=   UTR5_INTERNAL;
  longContentTypeMap["final-UTR5"]=      UTR5_FINAL;
  longContentTypeMap["single-UTR5"]=     UTR5_SINGLE;
  longContentTypeMap["initial-UTR3"]=    UTR3_INITIAL;
  longContentTypeMap["internal-UTR3"]=   UTR3_INTERNAL;
  longContentTypeMap["final-UTR3"]=      UTR3_FINAL;
  longContentTypeMap["single-UTR3"]=     UTR3_SINGLE;
}



Strand getStrand(ContentType t)
{
  switch(t)
    {
    case UNKNOWN_CONTENT_FORWARD:
    case INITIAL_EXON:
    case INTERNAL_EXON:
    case FINAL_EXON:
    case SINGLE_EXON:
    case INTRON:
    case UTR5_INTRON:
    case UTR3_INTRON:
    case INTERGENIC://###
    case UTR5_INITIAL:
    case UTR5_INTERNAL:
    case UTR5_FINAL:
    case UTR5_SINGLE:
    case UTR3_INITIAL:
    case UTR3_INTERNAL:
    case UTR3_FINAL:
    case UTR3_SINGLE:
      return FORWARD_STRAND;

    case UNKNOWN_CONTENT_REVERSE:
    case NEG_INITIAL_EXON:
    case NEG_INTERNAL_EXON:
    case NEG_FINAL_EXON:
    case NEG_SINGLE_EXON:
    case NEG_INTRON:
    case NEG_UTR5_INTRON:
    case NEG_UTR3_INTRON:
    case NEG_UTR5_INITIAL:
    case NEG_UTR5_INTERNAL:
    case NEG_UTR5_FINAL:
    case NEG_UTR5_SINGLE:
    case NEG_UTR3_INITIAL:
    case NEG_UTR3_INTERNAL:
    case NEG_UTR3_FINAL:
    case NEG_UTR3_SINGLE:
      return REVERSE_STRAND;

    default: throw BOOM::String(__FILE__)+__LINE__;
    }
}



BOOM::String contentTypeToString(ContentType t)
{
  switch(t)
    {
    case UNKNOWN_CONTENT_FORWARD: return "UNKNOWN-CONTENT-FORWARD";
    case UNKNOWN_CONTENT_REVERSE: return "UNKNOWN-CONTENT-REVERSE";
    case INITIAL_EXON:         return "INITIAL-EXON";
    case INTERNAL_EXON:        return "INTERNAL-EXON";
    case FINAL_EXON:           return "FINAL-EXON";
    case SINGLE_EXON:          return "SINGLE-EXON";
    case INTRON:               return "INTRON";
    case UTR5_INTRON:          return "UTR5-INTRON";
    case UTR3_INTRON:          return "UTR3-INTRON";
    case INTERGENIC:           return "INTERGENIC";

    case UTR5_INITIAL:         return "UTR5-INITIAL";
    case UTR5_INTERNAL:        return "UTR5-INTERNAL";
    case UTR5_FINAL:           return "UTR5-FINAL";
    case UTR5_SINGLE:          return "UTR5-SINGLE";
    case UTR3_INITIAL:         return "UTR3-INITIAL";
    case UTR3_INTERNAL:        return "UTR3-INTERNAL";
    case UTR3_FINAL:           return "UTR3-FINAL";
    case UTR3_SINGLE:          return "UTR3-SINGLE";

    case NEG_INITIAL_EXON:     return "NEG-INITIAL-EXON";
    case NEG_INTERNAL_EXON:    return "NEG-INTERNAL-EXON";
    case NEG_FINAL_EXON:       return "NEG-FINAL-EXON";
    case NEG_SINGLE_EXON:      return "NEG-SINGLE-EXON";
    case NEG_INTRON:           return "NEG-INTRON";
    case NEG_UTR5_INTRON:      return "NEG-UTR5-INTRON";
    case NEG_UTR3_INTRON:      return "NEG-UTR3-INTRON";

    case NEG_UTR5_INITIAL:         return "NEG-UTR5-INITIAL";
    case NEG_UTR5_INTERNAL:        return "NEG-UTR5-INTERNAL";
    case NEG_UTR5_FINAL:           return "NEG-UTR5-FINAL";
    case NEG_UTR5_SINGLE:          return "NEG-UTR5-SINGLE";
    case NEG_UTR3_INITIAL:         return "NEG-UTR3-INITIAL";
    case NEG_UTR3_INTERNAL:        return "NEG-UTR3-INTERNAL";
    case NEG_UTR3_FINAL:           return "NEG-UTR3-FINAL";
    case NEG_UTR3_SINGLE:          return "NEG-UTR3-SINGLE";

    default: throw BOOM::String("INVALID CONTENT TYPE: ")+int(t);
    }
}



BOOM::String contentTypeNiceString(ContentType t)
{
  switch(t)
    {
    case UNKNOWN_CONTENT_FORWARD: return "UNKNOWN-CONTENT-FORWARD";
    case UNKNOWN_CONTENT_REVERSE: return "UNKNOWN-CONTENT-REVERSE";
    case UTR5_INITIAL:
    case UTR3_INITIAL:
    case INITIAL_EXON:         return "initial-exon";
    case UTR5_INTERNAL:
    case UTR3_INTERNAL:
    case INTERNAL_EXON:        return "internal-exon";
    case UTR5_FINAL:
    case UTR3_FINAL:
    case FINAL_EXON:           return "final-exon";
    case UTR5_SINGLE:
    case UTR3_SINGLE:
    case SINGLE_EXON:          return "single-exon";
    case UTR5_INTRON: 
    case UTR3_INTRON:
    case INTRON:               return "intron";
    case INTERGENIC:           return "intergenic";
      //case UTR5_INTRON:          return "UTR5-intron";
      //case UTR3_INTRON:          return "UTR3-intron";
      //case UTR5_INITIAL:         return "initial-UTR5";
      //case UTR5_INTERNAL:        return "internal-UTR5";
      //case UTR5_FINAL:           return "final-UTR5";
      //case UTR5_SINGLE:          return "single-UTR5";
      //case UTR3_INITIAL:         return "initial-UTR3";
      //case UTR3_INTERNAL:        return "internal-UTR3";
      //case UTR3_FINAL:           return "final-UTR3";
      //case UTR3_SINGLE:          return "single-UTR3";
    case NEG_INITIAL_EXON:     return "initial-exon";
    case NEG_INTERNAL_EXON:    return "internal-exon";
    case NEG_FINAL_EXON:       return "final-exon";
    case NEG_SINGLE_EXON:      return "single-exon";
    case NEG_INTRON:           return "intron";
    case NEG_UTR5_INTRON:      return "UTR5-intron";
    case NEG_UTR3_INTRON:      return "UTR3-intron";
    case NEG_UTR5_INITIAL:         return "initial-UTR5";
    case NEG_UTR5_INTERNAL:        return "internal-UTR5";
    case NEG_UTR5_FINAL:           return "final-UTR5";
    case NEG_UTR5_SINGLE:          return "single-UTR5";
    case NEG_UTR3_INITIAL:         return "initial-UTR3";
    case NEG_UTR3_INTERNAL:        return "internal-UTR3";
    case NEG_UTR3_FINAL:           return "final-UTR3";
    case NEG_UTR3_SINGLE:          return "single-UTR3";

    default: throw BOOM::String("INVALID CONTENT TYPE: ")+int(t);
    }
}



ContentType stringToContentType(const BOOM::String &str)
{
  if(contentTypeMap.isDefined(str)) return contentTypeMap[str];
  if(longContentTypeMap.isDefined(str)) return longContentTypeMap[str];
  throw BOOM::String("undefined ContentType: ")+str;
}



istream &operator>>(istream &is,ContentType &t)
{
  BOOM::String buf;
  is >> buf;
  t=stringToContentType(buf);
  return is;
}



ostream &operator<<(ostream &os,const ContentType &t)
{
  os << contentTypeToString(t);
  return os;
}



ContentType reverseComplement(ContentType t)
{
  switch(t)
    {
    case INITIAL_EXON:         return NEG_INITIAL_EXON;
    case INTERNAL_EXON:        return NEG_INTERNAL_EXON;
    case FINAL_EXON:           return NEG_FINAL_EXON;
    case SINGLE_EXON:          return NEG_SINGLE_EXON;
    case INTRON:               return NEG_INTRON;
    case UTR5_INTRON:          return NEG_UTR5_INTRON;
    case UTR3_INTRON:          return NEG_UTR3_INTRON;
    case UTR5_INITIAL:         return NEG_UTR5_INITIAL;
    case UTR5_INTERNAL:        return NEG_UTR5_INTERNAL;
    case UTR5_FINAL:           return NEG_UTR5_FINAL;
    case UTR5_SINGLE:          return NEG_UTR5_SINGLE;
    case UTR3_INITIAL:         return NEG_UTR3_INITIAL;
    case UTR3_INTERNAL:        return NEG_UTR3_INTERNAL;
    case UTR3_FINAL:           return NEG_UTR3_FINAL;
    case UTR3_SINGLE:          return NEG_UTR3_SINGLE;
    case UNKNOWN_CONTENT_FORWARD: return UNKNOWN_CONTENT_REVERSE;

    case NEG_INITIAL_EXON:     return INITIAL_EXON;
    case NEG_INTERNAL_EXON:    return INTERNAL_EXON;
    case NEG_FINAL_EXON:       return FINAL_EXON;
    case NEG_SINGLE_EXON:      return SINGLE_EXON;
    case NEG_INTRON:           return INTRON;
    case NEG_UTR5_INTRON:      return UTR5_INTRON;
    case NEG_UTR3_INTRON:      return UTR3_INTRON;
    case NEG_UTR5_INITIAL:         return UTR5_INITIAL;
    case NEG_UTR5_INTERNAL:        return UTR5_INTERNAL;
    case NEG_UTR5_FINAL:           return UTR5_FINAL;
    case NEG_UTR5_SINGLE:          return UTR5_SINGLE;
    case NEG_UTR3_INITIAL:         return UTR3_INITIAL;
    case NEG_UTR3_INTERNAL:        return UTR3_INTERNAL;
    case NEG_UTR3_FINAL:           return UTR3_FINAL;
    case NEG_UTR3_SINGLE:          return UTR3_SINGLE;
    case UNKNOWN_CONTENT_REVERSE:  return UNKNOWN_CONTENT_FORWARD;
    case INTERGENIC:           return INTERGENIC;

    default: INTERNAL_ERROR;
    }
}



bool isCoding(ContentType t)
{
  switch(t)
    {
    case INITIAL_EXON:         
    case INTERNAL_EXON:        
    case FINAL_EXON:           
    case SINGLE_EXON:          
    case NEG_INITIAL_EXON:     
    case NEG_INTERNAL_EXON:    
    case NEG_FINAL_EXON:       
    case NEG_SINGLE_EXON:      
      return true;

    case INTRON:  
    case UTR5_INTRON:
    case UTR3_INTRON:
    case INTERGENIC:           
    case NEG_INTRON:  
    case NEG_UTR5_INTRON:
    case NEG_UTR3_INTRON:
    case UTR5_INITIAL:
    case UTR5_INTERNAL:
    case UTR5_FINAL:
    case UTR5_SINGLE:
    case UTR3_INITIAL:
    case UTR3_INTERNAL:
    case UTR3_FINAL:
    case UTR3_SINGLE:
    case NEG_UTR5_INITIAL:
    case NEG_UTR5_INTERNAL:
    case NEG_UTR5_FINAL:
    case NEG_UTR5_SINGLE:
    case NEG_UTR3_INITIAL:
    case NEG_UTR3_INTERNAL:
    case NEG_UTR3_FINAL:
    case NEG_UTR3_SINGLE:
      return false;

    case UNKNOWN_CONTENT_FORWARD:
    case UNKNOWN_CONTENT_REVERSE:
      throw "isCoding(UNKNOWN_CONTENT_FORWARD/REVERSE)";

    default: INTERNAL_ERROR;
    }
}



bool isIntron(ContentType t)
{
  switch(t)
    {
    case INTRON:
    case UTR5_INTRON:
    case UTR3_INTRON:
    case NEG_INTRON: 
    case NEG_UTR5_INTRON:
    case NEG_UTR3_INTRON:
      return true;
    default:
      return false;
    }
}



bool isIntergenic(ContentType t)
{
  return t==INTERGENIC;
}



bool isUTR(ContentType t)
{
  switch(t)
    {
    case UTR5_INITIAL:
    case UTR5_INTERNAL:
    case UTR5_FINAL:
    case UTR5_SINGLE:
    case UTR3_INITIAL:
    case UTR3_INTERNAL:
    case UTR3_FINAL:
    case UTR3_SINGLE:
    case NEG_UTR5_INITIAL:
    case NEG_UTR5_INTERNAL:
    case NEG_UTR5_FINAL:
    case NEG_UTR5_SINGLE:
    case NEG_UTR3_INITIAL:
    case NEG_UTR3_INTERNAL:
    case NEG_UTR3_FINAL:
    case NEG_UTR3_SINGLE:
      return true;

    default: 
      return false;
    }
}



bool isUTR5(ContentType t)
{
  switch(t)
    {
    case UTR5_INITIAL:
    case UTR5_INTERNAL:
    case UTR5_FINAL:
    case UTR5_SINGLE:
    case NEG_UTR5_INITIAL:
    case NEG_UTR5_INTERNAL:
    case NEG_UTR5_FINAL:
    case NEG_UTR5_SINGLE:
      return true;
    }
  return false;
}



bool isUTR3(ContentType t)
{
  switch(t)
    {
    case UTR3_INITIAL:
    case UTR3_INTERNAL:
    case UTR3_FINAL:
    case UTR3_SINGLE:
    case NEG_UTR3_INITIAL:
    case NEG_UTR3_INTERNAL:
    case NEG_UTR3_FINAL:
    case NEG_UTR3_SINGLE:
      return true;
    }
}



bool isUTR5intron(ContentType t)
{
  switch(t)
    {
    case UTR5_INTRON:
    case NEG_UTR5_INTRON:
      return true;
    }
}



bool isUTR3intron(ContentType t)
{
  switch(t)
    {
    case UTR3_INTRON:
    case NEG_UTR3_INTRON:
      return true;
    }
}



static Set<SignalType> sigset(SignalType t)
{
  Set<SignalType> s;
  s+=t;
  return s;
}


static Set<SignalType> sigset(SignalType t1,SignalType t2)
{
  Set<SignalType> s;
  s+=t1;
  s+=t2;
  return s;
}


static Set<SignalType> sigset(SignalType t1,SignalType t2,SignalType t3)
{
  Set<SignalType> s;
  s+=t1;
  s+=t2;
  s+=t3;
  return s;
}


Set<SignalType> leftSignals(ContentType contentType)
{
  switch(contentType)
    {
    case INTERGENIC:               return sigset(LEFT_TERMINUS,TAG,TES);
    case INITIAL_EXON:             return sigset(ATG);
    case INTERNAL_EXON:            return sigset(AG);
    case FINAL_EXON:               return sigset(AG);
    case SINGLE_EXON:              return sigset(ATG);
    case INTRON:                   return sigset(GT);
    case UTR5_INTRON:              return sigset(UTR5GT);
    case UTR3_INTRON:              return sigset(UTR3GT);
    case UTR5_INITIAL:             return sigset(TSS);
    case UTR5_INTERNAL:            return sigset(UTR5AG);
    case UTR5_FINAL:               return sigset(UTR5AG);
    case UTR5_SINGLE:              return sigset(TSS);
    case UTR3_INITIAL:             return sigset(TAG);
    case UTR3_INTERNAL:            return sigset(UTR3AG);
    case UTR3_FINAL:               return sigset(UTR3AG);
    case UTR3_SINGLE:              return sigset(TAG);
    case NEG_INITIAL_EXON:         return sigset(NEG_GT);
    case NEG_INTERNAL_EXON:        return sigset(NEG_GT);
    case NEG_FINAL_EXON:           return sigset(NEG_TAG);
    case NEG_SINGLE_EXON:          return sigset(NEG_TAG);
    case NEG_INTRON:               return sigset(NEG_AG);
    case NEG_UTR5_INTRON:          return sigset(NEG_UTR5AG);
    case NEG_UTR3_INTRON:          return sigset(NEG_UTR3AG);
    case NEG_UTR5_INITIAL:         return sigset(NEG_GT);
    case NEG_UTR5_INTERNAL:        return sigset(NEG_GT);
    case NEG_UTR5_FINAL:           return sigset(NEG_ATG,NEG_TES);
    case NEG_UTR5_SINGLE:          return sigset(NEG_ATG,NEG_TES);
    case NEG_UTR3_INITIAL:         return sigset(NEG_GT);
    case NEG_UTR3_INTERNAL:        return sigset(NEG_GT);
    case NEG_UTR3_FINAL:           return sigset(NEG_TES);
    case NEG_UTR3_SINGLE:          return sigset(NEG_TES);
    default: INTERNAL_ERROR;
    }
}



Set<SignalType> rightSignals(ContentType contentType)
{
  switch(contentType)
    {
    case INTERGENIC:               return sigset(TSS,ATG,RIGHT_TERMINUS);
    case INITIAL_EXON:             return sigset(GT);
    case INTERNAL_EXON:            return sigset(GT);
    case FINAL_EXON:               return sigset(TAG);
    case SINGLE_EXON:              return sigset(TAG);
    case INTRON:                   return sigset(AG);
    case UTR5_INTRON:              return sigset(UTR5AG);
    case UTR3_INTRON:              return sigset(UTR3AG);
    case UTR5_INITIAL:             return sigset(UTR5GT);
    case UTR5_INTERNAL:            return sigset(UTR5GT);
    case UTR5_FINAL:               return sigset(ATG,TES);
    case UTR5_SINGLE:              return sigset(ATG,TES);
    case UTR3_INITIAL:             return sigset(UTR3GT);
    case UTR3_INTERNAL:            return sigset(UTR3GT);
    case UTR3_FINAL:               return sigset(TES);
    case UTR3_SINGLE:              return sigset(TES);
    case NEG_INITIAL_EXON:         return sigset(NEG_ATG);
    case NEG_INTERNAL_EXON:        return sigset(NEG_AG);
    case NEG_FINAL_EXON:           return sigset(NEG_AG);
    case NEG_SINGLE_EXON:          return sigset(NEG_ATG);
    case NEG_INTRON:               return sigset(NEG_GT);
    case NEG_UTR5_INTRON:          return sigset(NEG_UTR5GT);
    case NEG_UTR3_INTRON:          return sigset(NEG_UTR3GT);
    case NEG_UTR5_INITIAL:         return sigset(NEG_TSS);
    case NEG_UTR5_INTERNAL:        return sigset(NEG_UTR5AG);
    case NEG_UTR5_FINAL:           return sigset(NEG_UTR5AG);
    case NEG_UTR5_SINGLE:          return sigset(NEG_TSS);
    case NEG_UTR3_INITIAL:         return sigset(NEG_TAG);
    case NEG_UTR3_INTERNAL:        return sigset(NEG_UTR3AG);
    case NEG_UTR3_FINAL:           return sigset(NEG_UTR3AG);
    case NEG_UTR3_SINGLE:          return sigset(NEG_TAG);
    default: throw BOOM::String(__FILE__)+__LINE__;
    }
}



ExonType contentTypeToExonType(ContentType ct)
{
  switch(ct) {
    case INITIAL_EXON:             return ET_INITIAL_EXON;
    case INTERNAL_EXON:            return ET_INTERNAL_EXON;
    case FINAL_EXON:               return ET_FINAL_EXON;
    case SINGLE_EXON:              return ET_SINGLE_EXON;
    case NEG_INITIAL_EXON:         return ET_INITIAL_EXON;
    case NEG_INTERNAL_EXON:        return ET_INTERNAL_EXON;
    case NEG_FINAL_EXON:           return ET_FINAL_EXON;
    case NEG_SINGLE_EXON:          return ET_SINGLE_EXON;
    case UTR5_INITIAL:             return ET_INITIAL_UTR5;
    case UTR5_INTERNAL:            return ET_INTERNAL_UTR5;
    case UTR5_FINAL:               return ET_FINAL_UTR5;
    case UTR5_SINGLE:              return ET_SINGLE_UTR5;
    case UTR3_INITIAL:             return ET_INITIAL_UTR3;
    case UTR3_INTERNAL:            return ET_INTERNAL_UTR3;
    case UTR3_FINAL:               return ET_FINAL_UTR3;
    case UTR3_SINGLE:              return ET_SINGLE_UTR3;
    case NEG_UTR5_INITIAL:         return ET_INITIAL_UTR5;
    case NEG_UTR5_INTERNAL:        return ET_INTERNAL_UTR5;
    case NEG_UTR5_FINAL:           return ET_FINAL_UTR5;
    case NEG_UTR5_SINGLE:          return ET_SINGLE_UTR5;
    case NEG_UTR3_INITIAL:         return ET_INITIAL_UTR3;
    case NEG_UTR3_INTERNAL:        return ET_INTERNAL_UTR3;
    case NEG_UTR3_FINAL:           return ET_FINAL_UTR3;
    case NEG_UTR3_SINGLE:          return ET_SINGLE_UTR3;
  default: INTERNAL_ERROR;
  }
}



ContentType exonTypeToContentType(ExonType exonType,
				  Strand strand)
{
  switch(exonType)
    {
    case ET_INITIAL_EXON:
      return strand==FORWARD_STRAND ? INITIAL_EXON : NEG_INITIAL_EXON;
    case ET_INTERNAL_EXON:
      return strand==FORWARD_STRAND ? INTERNAL_EXON : NEG_INTERNAL_EXON;
    case ET_FINAL_EXON:
      return strand==FORWARD_STRAND ? FINAL_EXON : NEG_FINAL_EXON;
    case ET_SINGLE_EXON:
      return strand==FORWARD_STRAND ? SINGLE_EXON : NEG_SINGLE_EXON;
    case ET_INITIAL_UTR5:
      return strand==FORWARD_STRAND ? UTR5_INITIAL : NEG_UTR5_INITIAL;
    case ET_INTERNAL_UTR5:
      return strand==FORWARD_STRAND ? UTR5_INTERNAL : NEG_UTR5_INTERNAL;
    case ET_FINAL_UTR5:
      return strand==FORWARD_STRAND ? UTR5_FINAL : NEG_UTR5_FINAL;
    case ET_SINGLE_UTR5:
      return strand==FORWARD_STRAND ? UTR5_SINGLE : NEG_UTR5_SINGLE;
    case ET_INITIAL_UTR3:
      return strand==FORWARD_STRAND ? UTR3_INITIAL : NEG_UTR3_INITIAL;
    case ET_INTERNAL_UTR3:
      return strand==FORWARD_STRAND ? UTR3_INTERNAL : NEG_UTR3_INTERNAL;
    case ET_FINAL_UTR3:
      return strand==FORWARD_STRAND ? UTR3_FINAL : NEG_UTR3_FINAL;
    case ET_SINGLE_UTR3:
      return strand==FORWARD_STRAND ? UTR3_SINGLE : NEG_UTR3_SINGLE;

    case ET_UTR: // nonspecific
    case ET_UTR5:// nonspecific
    case ET_UTR3:// nonspecific
    case ET_EXON:// nonspecific
    default:
      throw "::exonTypeToContentType()";
    }
}



ContentType dropStrand(ContentType t)
{
  switch(t)
    {
    case INITIAL_EXON:
    case INTERNAL_EXON:
    case FINAL_EXON:
    case SINGLE_EXON:
    case INTRON:
    case UTR5_INTRON:
    case UTR3_INTRON:
    case UTR5_INITIAL:
    case UTR5_INTERNAL:
    case UTR5_FINAL:
    case UTR5_SINGLE:
    case UTR3_INITIAL:
    case UTR3_INTERNAL:
    case UTR3_FINAL:
    case UTR3_SINGLE:
    case INTERGENIC:           
      return t;
    case UNKNOWN_CONTENT_FORWARD:
    case NEG_INITIAL_EXON:
    case NEG_INTERNAL_EXON:
    case NEG_FINAL_EXON:
    case NEG_SINGLE_EXON:
    case NEG_INTRON:
    case NEG_UTR5_INTRON:
    case NEG_UTR3_INTRON:
    case NEG_UTR5_INITIAL:
    case NEG_UTR5_INTERNAL:
    case NEG_UTR5_FINAL:
    case NEG_UTR5_SINGLE:
    case NEG_UTR3_INITIAL:
    case NEG_UTR3_INTERNAL:
    case NEG_UTR3_FINAL:
    case NEG_UTR3_SINGLE:
    case UNKNOWN_CONTENT_REVERSE:
      return reverseComplement(t);
    }
  INTERNAL_ERROR;
}


