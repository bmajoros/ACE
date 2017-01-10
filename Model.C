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
    spliceShiftDistr(NULL)
{
  // ctor
}



Model::~Model()
{
  delete exonLengthDistr;
  delete intronLengthDistr;
  delete intergenicLengthDistr;
  delete spliceShiftDistr;

  /* These are allocated on the stack -- DO NOT DELETE THEM!
  delete signalSensors;
  delete contentSensors;
  */
}




