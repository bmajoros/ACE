/****************************************************************
 ACEplus.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ACEplus_H
#define INCL_ACEplus_H
#include "ACE.H"
#include "ContentSensors.H"
#include "Model.H"
#include "TranscriptPaths.H"
using namespace std;
using namespace BOOM;


class ACEplus : public ACE {
public:
  ACEplus();
  virtual ~ACEplus() {}
  virtual int main(int argc,char *argv[]);
protected:
  ContentSensors contentSensors;
  Model model;
  virtual void parseCommandLine(const CommandLine &);
  virtual void processConfig(const String &filename);
  virtual ContentSensor *loadContentSensor(const String &label,
					   ConfigFile &);
  virtual void buildPSAs(ContentSensors &contentSensors,int seqLen,
			 Sequence &seq,String &str);
  virtual void buildPSA(ContentType type,ContentSensors &contentSensors,
			int seqLen,Sequence &seq,String &str);
  void checkProjection(const String &outGff,
		       bool &mapped,
		       const Labeling &refLab,
		       const Labeling &projectedLab,
		       ostream &osACE);
  virtual void enumerateAlts(TranscriptPaths &,
			     Essex::CompositeNode *altTransEssex,
			     TranscriptSignals *signals,
			     GffTranscript *altTrans,
			     ostream &osACE,
			     const Labeling &refLab,
			     const Labeling &projectedLab);
  void processAltStructure(TranscriptPath &,
			   Essex::CompositeNode *altStructNode,
			   const Labeling &projectedLab,
			   int whichStructure,
			   TranscriptSignals &);
  virtual Essex::CompositeNode * processAltStructure(AlternativeStructure &,
				   Essex::CompositeNode *altStructNode,
				   const Labeling &projectedLab);
  void refineStartCodon(int start,GffTranscript &,Essex::CompositeNode *&msg);
  double getRefLikelihood(const Labeling &refLab,GffTranscript *altTrans);
};

#endif

