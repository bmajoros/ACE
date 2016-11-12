/****************************************************************
 IsochoreTable.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_IsochoreTable_H
#define INCL_IsochoreTable_H
#include "BOOM/Regex.H"
#include "BOOM/Vector.H"
#include "BOOM/String.H"
#include "BOOM/ConfigFile.H"
#include "BOOM/Map.H"
#include "DiscreteDistribution.H"
#include "ContentType.H"
#include "Isochore.H"
#include "GarbageCollector.H"

class IsochoreTable {
  static BOOM::Regex lineRegex, globalRegex, commentRegex;
  BOOM::Vector<Isochore*> isochores;
  BOOM::Vector<float> boundaries;
  int queueCapacity;
  GarbageCollector &garbageCollector;
  static BOOM::String romanNumerals[21];

  void loadIsochore(BOOM::ConfigFile &,Isochore &);
  void loadSubmodels(BOOM::ConfigFile &configFile,Isochore &);
  void loadContentSensor(BOOM::ConfigFile &configFile,
			 const BOOM::String &modelLabel,
			 DistributionType distrType,
			 const BOOM::String &distrLabel,
			 ContentType contentType,
			 Isochore &);
  void loadSignalSensor(BOOM::ConfigFile &configFile,
			const BOOM::String &modelLabel,
			const BOOM::String &consensusLabel,
			Isochore &);
  void createVirtualSensor(SignalType,Isochore &);
  void loadSignalPeptideModel(const BOOM::String &filename,
			      BOOM::ConfigFile &configFile,Isochore &);
  BOOM::String toRomanNumerals(int);
public:
  IsochoreTable(GarbageCollector &);
  void load(const BOOM::String &filename);
  Isochore *getIsochore(float GC);
  Isochore *getIthIsochore(int i);
  int getQueueCapacity();
  int getNumIsochores();
};


#endif
