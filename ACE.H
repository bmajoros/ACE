/****************************************************************
 ACE.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ACE_H
#define INCL_ACE_H
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/CigarString.H"
#include "BOOM/GffReader.H"
#include "BOOM/FastaReader.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/CombinationIterator.H"
#include "BOOM/Essex.H"
#include "BOOM/Regex.H"
#include "BOOM/ConfigFile.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/BandedSmithWaterman.H"
#include "Labeling.H"
#include "ProjectionChecker.H"
#include "EnumerateAltStructures.H"
#include "SignalSensors.H"
#include "GarbageCollector.H"
#include "StartCodonFinder.H"
#include "NMD.H"
#include "Variant.H"
#include "VariantClassifier.H"
#include "OrfAnalyzer.H"
using namespace std;
using namespace BOOM;

/****************************************************************
 class ACE
 ****************************************************************/
class ACE {
public:
  ACE();
  virtual ~ACE();
  virtual int main(int argc,char *argv[]);
protected:
  NMD nmd;
  FastaWriter fastaWriter;
  SignalSensors sensors;
  String labelingFile;
  String substrate, altDefline, xmlFilename, globalCoord;
  int MAX_SPLICE_SHIFT, MIN_EXON_LEN, MIN_INTRON_LEN, NMD_DISTANCE_PARM;
  bool allowExonSkipping, allowIntronRetention, allowCrypticSites;
  bool reverseStrand, quiet;
  String CIGAR;
  CigarString cigar;
  OrfAnalyzer *orfAnalyzer;
  Essex::CompositeNode *root, *status;
  Essex::Node *startCodonMsg;
  Vector<Variant> variants;
  Regex warningsRegex, errorsRegex, variantRegex, coordRegex;
  int VCFwarnings, VCFerrors;
  String refSeqStr, altSeqStr;
  Sequence refSeq, altSeq;
  int refSeqLen, altSeqLen;
  GffTranscript *refTrans;
  String refProtein;
  double openPenalty, extendPenalty;
  int bandwidth;
  CigarAlignment *alignment, *revAlignment;
  SubstitutionMatrix<float> *substMatrix; // protein matrix
  GarbageIgnorer garbageCollector;
  String configFile, refGffFile, refFasta, altFasta, outGff, outACE;

  virtual void parseCommandLine(const CommandLine &);
  virtual void commandLineOpts(const CommandLine &);
  virtual void loadInputs(const String &configFile,
		  const String &refGffFile,
		  const String &refFasta,
		  const String &altFasta);
  virtual void buildAlignment();
  virtual GffTranscript *loadGff(const String &filename);
  virtual String loadSeq(const String &filename);
  virtual String loadSeq(const String &filename,String &cigar);
  virtual void computeLabeling(GffTranscript &,Labeling &);
  virtual void mapLabeling(Labeling &from,Labeling &to,const String &filename);
  virtual void mapTranscript(const String &outfile);
  virtual bool mapExon(GffExon &,CigarAlignment &);
  virtual void writeProtein(const String &defline,
		    const String &protein,
		    const String &filename);
  virtual void append(Essex::CompositeNode *,const String &tag,
		      const String &message);
  virtual void append(Essex::CompositeNode *,const char *tag,
		      const char *message);
  virtual void append(Essex::CompositeNode *,const char *tag,
		      const String &message);
  virtual void append(Essex::CompositeNode *,const char *tag,int);
  virtual void appendBrokenSignals(const TranscriptSignals *);
  virtual void writeXML();
  virtual void processConfig(const String &filename);
  virtual void parseConsensusList(const String &tag,
				  ConfigFile &,
				  Set<String> &into);
  virtual SignalSensor *loadModel(const String &label,
				  ConfigFile &);
  virtual double alignProteins(const String &refStr,
			      const String &altStr,
			      int &matches);
  virtual void percentMatch(int matches,
			    int refLen,
			    int altLen,
			    Essex::CompositeNode *parent);
  virtual void parseVariants(const String &,
			     Vector<Variant> &,
			     int substrateLen);
  virtual Essex::CompositeNode *makeEssexVariants();
  virtual bool checkRefGene();
  virtual void checkProjection(const String &outGff,
			       bool &mapped,
			       const Labeling &projectedLab,
			       ostream &osACE);
  virtual void handleProteinFate(const AlternativeStructure &,
				 Essex::CompositeNode *);
  virtual void enumerateAlts(Essex::CompositeNode *altTransEssex,
			     TranscriptSignals *signals,
			     GffTranscript *altTrans,
			     ostream &osACE,
			     const Labeling &projectedLab);
  virtual Essex::CompositeNode *processAltStructure(AlternativeStructure &,
				   Essex::CompositeNode *altStructNode,
				   const Labeling &projectedLab);
  virtual void listStructureChanges(const AlternativeStructure &,
				    Essex::CompositeNode *,
				    Essex::CompositeNode *&msg,
				    Essex::CompositeNode *gainOfCoding);
  virtual void handleNoncoding(const GffTranscript *altTrans);
  virtual void handleCoding(GffTranscript *altTrans,
			    ProjectionChecker &checker,
			    const Labeling &projectedLab);
  virtual void initEssex(ostream &osACE,
			 const CommandLine &);
  virtual void flushOutput(ostream &osACE,const bool &mapped);
  virtual int getTruncationLength(const GffTranscript &transcript,
				  int PTC,
				  int stop);
  virtual Essex::Node *signalNode(const String &tag,
				  const String &signal,
				  double score,
				  double cutoff);
  virtual void analyzeEarlierStart(GffTranscript *altTrans,
				   ProjectionChecker &,
				   const Labeling &projectedLab,
				   Essex::CompositeNode *status);
};


#endif
