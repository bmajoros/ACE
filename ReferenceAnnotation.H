/****************************************************************
 ReferenceAnnotation.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ReferenceAnnotation_H
#define INCL_ReferenceAnnotation_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/Comparator.H"
#include "Labeling.H"
#include "LabelMatrix.H"
#include "Signal.H"
#include "ContentRegions.H"
#include "Isochore.H"

using namespace std;
using namespace BOOM;


class ReferenceAnnotation {
public:
  ReferenceAnnotation();
  ReferenceAnnotation(const String &annoGFFfile,const String &labelingFile,
		      const String &matrixFile,Isochore &,
		      const String &altSeqStr,const Sequence &altSequence);
  virtual ~ReferenceAnnotation();
  const LabelMatrix &getMatrix() const;
  const Labeling &getLabeling() const;
  const ContentRegions &getRegions() const;
  const Vector<Signal*> &getSignals() const;
  SignalPtr getStartCodon() const;
  int getStartPosition() const;
  const String &getAltSeqStr() { return altSeqStr; }
  const Sequence &getAltSeq() { return altSeq; }
protected:
  const String altSeqStr;
  const Sequence altSeq;
  LabelMatrix *matrix;
  Labeling labeling;
  Vector<Signal*> signals;
  ContentRegions *contentRegions;
  void loadMatrix(const String &filename,const String &baseDir);
  void loadLabeling(const String &filename);
  void loadGFF(const String &filename,const int seqLen);
  void initSignals(Isochore &,const String &altSeqStr,const Sequence &altSequence);
  void makeSignal(SignalType,int intervalBeginOrEnd,const String &altSequence,
		  const Sequence &,SignalSensor &);
  void makeSignal(SignalType,const String &,const Sequence &,int contextWindowPos,
		  SignalSensor &);
  void sortSignals();
};

#endif

