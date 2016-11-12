/****************************************************************
 TopologyLoader.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_TopologyLoader_H
#define INCL_TopologyLoader_H
#include "BOOM/Regex.H"
#include "SignalTypeProperties.H"

class TopologyLoader
{
  static BOOM::Regex transitionRegex;      // "GT -> AG : INTRON"
  static BOOM::Regex phaseRegex;           // "GT has phase 0"
  static BOOM::Regex strandRegex;          // "GT strand +"
  static BOOM::Regex consensusCodingRegex; // "ATG consensus coding"
  static BOOM::Regex commentRegex;         // "# comment"
  static BOOM::Regex enableRegex;          // "enable ATG"
  static BOOM::Regex blankLineRegex;
public:
  static SignalTypeProperties &load(const BOOM::String &filename);
};


#endif
