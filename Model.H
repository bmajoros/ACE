/****************************************************************
 Model.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Model_H
#define INCL_Model_H
#include <iostream>
#include "SignalSensors.H"
#include "ContentSensors.H"
#include "DiscreteDistribution.H"
using namespace std;
using namespace BOOM;

struct Model {
  SignalSensors *signalSensors;
  ContentSensors *contentSensors;
  DiscreteDistribution *exonLengthDistr;
  DiscreteDistribution *intronLengthDistr;
  DiscreteDistribution *intergenicLengthDistr;
  DiscreteDistribution *spliceShiftDistr;
  int MAX_SPLICE_SHIFT, MIN_EXON_LEN, MIN_INTRON_LEN, NMD_DISTANCE_PARM;
  int MAX_CRYPTIC_EXON_LEN; // maximum length of a cryptic exon
  float EXON_STRENGTHENING_THRESHOLD; // min ratio of alt/ref for intronic
                                      // sequence to (possibly) become exonic
  float EXON_WEAKENING_THRESHOLD; // below this alt/ref ratio, exonic sequence
                                  // can become intronic
  float MIN_EXON_INTRON_RATIO; // ratio of scores for exon/intron sensors;
                               // below this, exon definition fails
  bool allowExonSkipping, allowIntronRetention, allowCrypticSites;
  bool allowDeNovoSites, allowCrypticExons, allowRegulatoryChanges;
  Model();
  virtual ~Model();
};

#endif

