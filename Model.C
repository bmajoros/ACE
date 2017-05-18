/****************************************************************
 Model.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Model.H"
using namespace std;
using namespace BOOM;


Model::Model()
  : signalSensors(NULL), contentSensors(NULL), exonLengthDistr(NULL),
    intronLengthDistr(NULL), intergenicLengthDistr(NULL),
    spliceShiftDistr(NULL), transitions(NULL),
    maxIntronRetentionLen(0), minIntronRetentionLLR(0.0),
    maxDeNovoExonLen(0), minDeNovoExonLLR(0.0), sensorScale(1.0)
{
  // ctor
}



Model::~Model()
{
  delete transitions;
  delete exonLengthDistr;
  delete intronLengthDistr;
  delete intergenicLengthDistr;
  delete spliceShiftDistr;

  /* These are allocated on the stack -- DO NOT DELETE THEM!
  delete signalSensors;
  delete contentSensors;
  */
}




