/****************************************************************
 TopologyLoader.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "TopologyLoader.H"
#include <iostream>
#include "BOOM/File.H"


// "GT -> AG : INTRON"
BOOM::Regex TopologyLoader::transitionRegex(
   "^\\s*(\\S+)\\s*->\\s*(\\S+)\\s*:\\s*(\\S+)\\s*$");

// "GT has phase 0"
BOOM::Regex TopologyLoader::phaseRegex(
   "^\\s*(\\S+)\\s+has\\s+phase\\s+(\\S+)\\s*$");

// "GT strand +"
BOOM::Regex TopologyLoader::strandRegex(
   "^\\s*(\\S+)\\s+strand\\s+(\\S+)\\s*$");

// "ATG consensus coding"
BOOM::Regex TopologyLoader::consensusCodingRegex(
   "^\\s*(\\S+)\\s+consensus\\s+coding\\s*$");

// "enable ATG"
BOOM::Regex TopologyLoader::enableRegex(
   "^\\s*enable\\s+(\\S+)$");

// "# comment"
BOOM::Regex TopologyLoader::commentRegex("^\\s*#.*$");

BOOM::Regex TopologyLoader::blankLineRegex("^\\s*$");



SignalTypeProperties &TopologyLoader::load(const BOOM::String &filename)
{
  SignalTypeProperties &stp=SignalTypeProperties::global;
  int lineNumber=0;
  BOOM::File file(filename);
  while(!file.eof()) {
    BOOM::String line=file.readLine();
    ++lineNumber;
    if(commentRegex.match(line)) continue;
    if(transitionRegex.search(line)) {
      // "GT -> AG : INTRON"
      SignalType from=stringToSignalType(transitionRegex[1]);
      SignalType to=stringToSignalType(transitionRegex[2]);
      ContentType content=stringToContentType(transitionRegex[3]);
      stp.belongTo(from,content);
      stp.linkTo(to,content);
    }
    else if(phaseRegex.search(line)) {
      // "GT has phase 0"
      SignalType sig=stringToSignalType(phaseRegex[1]);
      int phase=phaseRegex[2].asInt();
      stp.allowPhase(sig,phase);
    }
    else if(strandRegex.search(line)) {
      // "GT strand +"
      SignalType sig=stringToSignalType(strandRegex[1]);
      char strand=strandRegex[2][0];
      stp.setStrand(sig,strand);
    }
    else if(consensusCodingRegex.search(line)){
      // "ATG consensus coding"
      SignalType sig=stringToSignalType(consensusCodingRegex[1]);
      stp.setConsensusCoding(sig,true);
    }
    else if(enableRegex.search(line)){
      // "enable ATG"
      SignalType sig=stringToSignalType(enableRegex[1]);
      stp.enable(sig);
    }
    else if(!commentRegex.match(line) && !blankLineRegex.match(line))
      throw filename+", line "+BOOM::String(lineNumber)+", syntax error:\n"+line;
  }
  stp.belongTo(RIGHT_TERMINUS,INTERGENIC);
  
  return stp;
}


