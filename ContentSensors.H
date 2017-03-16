/****************************************************************
 ContentSensors.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ContentSensors_H
#define INCL_ContentSensors_H
#include <iostream>
#include "BOOM/Set.H"
#include "BOOM/String.H"
#include "ContentSensor.H"
#include "PrefixSumArray.H"
using namespace std;
using namespace BOOM;


class ContentSensors {
public:
  ContentSensors();
  virtual ~ContentSensors();
  void setSensor(ContentType,ContentSensor *);
  void setSpliceBackgroundModel(ContentSensor *);
  ContentSensor *getSpliceBackground();
  const ContentSensor *getSensor(ContentType) const;
  PrefixSumArray &getPSA(ContentType);
  double score(ContentType,int begin,int end) const;
protected:
  ContentSensor *exonSensor;
  ContentSensor *intronSensor;
  ContentSensor *intergenicSensor;
  ContentSensor *spliceBackgroundModel;
  PrefixSumArray exonPSA, intronPSA, intergenicPSA;
  ContentType simplifyContentType(ContentType type);
};

#endif

