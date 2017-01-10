/****************************************************************
 ContentSensors.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ContentSensors.H"
using namespace std;
using namespace BOOM;



ContentSensors::ContentSensors()
  : exonSensor(NULL), intronSensor(NULL), intergenicSensor(NULL)
{
  //ctor
}



const ContentSensor *ContentSensors::getSensor(ContentType type) const
{
  switch(simplifyContentType(type)) {
  case EXON:
    return exonSensor;
  case INTRON:
    return intronSensor;
  case INTERGENIC:
    return intergenicSensor;
  }
  INTERNAL_ERROR;
}



PrefixSumArray &ContentSensors::getPSA(ContentType type)
{
  switch(simplifyContentType(type)) {
  case EXON:
    return exonPSA;
  case INTRON:
    return intronPSA;
  case INTERGENIC:
    return intergenicPSA;
  }
  INTERNAL_ERROR;
}



double ContentSensors::score(ContentType type,int begin,int end) const
{
  const PrefixSumArray &psa=getPSA(type);
  return psa.getInterval(begin,end);
}



void ContentSensors::setSensor(ContentType type,ContentSensor *sensor)
{
  sensor->setContentType(type);
  switch(simplifyContentType(type)) {
  case SINGLE_EXON: // fall through...
  case EXON:           exonSensor=sensor;       return;
  case INTRON:         intronSensor=sensor;     return;
  case INTERGENIC:     intergenicSensor=sensor; return;
  }
  cout<<simplifyContentType(type)<<endl;
  INTERNAL_ERROR;
}



ContentType ContentSensors::simplifyContentType(ContentType type)
{
  switch(type) {
  case INITIAL_EXON:
  case INTERNAL_EXON:
  case FINAL_EXON:
  case SINGLE_EXON:
  case UTR5_INITIAL:
  case UTR5_INTERNAL:
  case UTR5_FINAL:
  case UTR5_SINGLE: // also serves as general 5' UTR marker
  case UTR3_INITIAL:
  case UTR3_INTERNAL:
  case UTR3_FINAL:
  case UTR3_SINGLE: // also serves as general 3' UTR marker
  case NEG_INITIAL_EXON:
  case NEG_INTERNAL_EXON:
  case NEG_FINAL_EXON:
  case NEG_SINGLE_EXON:
  case NEG_UTR5_INITIAL:
  case NEG_UTR5_INTERNAL:
  case NEG_UTR5_FINAL:
  case NEG_UTR5_SINGLE:
  case NEG_UTR3_INITIAL:
  case NEG_UTR3_INTERNAL:
  case NEG_UTR3_FINAL:
  case NEG_UTR3_SINGLE:
    return EXON;

  case INTRON:
  case UTR5_INTRON:
  case UTR3_INTRON:
  case NEG_INTRON:
  case NEG_UTR5_INTRON:
  case NEG_UTR3_INTRON:
    return INTRON;

  case INTERGENIC:
    return INTERGENIC;
  }
  INTERNAL_ERROR;
}

