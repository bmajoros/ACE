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
#include "Transitions.H"
using namespace std;
using namespace BOOM;

struct Model {
  Transitions *transitions;
  SignalSensors *signalSensors;
  ContentSensors *contentSensors;
  DiscreteDistribution *exonLengthDistr;
  DiscreteDistribution *intronLengthDistr;
  DiscreteDistribution *intergenicLengthDistr;
  DiscreteDistribution *spliceShiftDistr;
  int maxIntronRetentionLen;
  float minIntronRetentionLLR;
  int maxDeNovoExonLen;
  float minDeNovoExonLLR;
  int MAX_SPLICE_SHIFT, MIN_EXON_LEN, MIN_INTRON_LEN, NMD_DISTANCE_PARM;
  float EXON_STRENGTHENING_THRESHOLD; // min ratio of alt/ref for intronic
                                      // sequence to (possibly) become exonic
  float EXON_WEAKENING_THRESHOLD; // below this alt/ref ratio, exonic sequence
                                  // can become intronic
  float MIN_EXON_INTRON_RATIO; // ratio of scores for exon/intron sensors;
                               // below this, exon definition fails
  float MIN_SCORE; // don't report any structure scoring less than this
  int MAX_ALT_STRUCTURES; // maximum number of alternative structures to
                          // predict when there are splicing changes
  float coefDenovoExon;   // adjustment to transition prob for de novo exons
  bool allowExonSkipping, allowIntronRetention, allowCrypticSites;
  bool allowDeNovoSites, allowCrypticExons, allowRegulatoryChanges;
  float sensorScale;
  float exonIntercept;
  Model();
  virtual ~Model();
};

#endif

