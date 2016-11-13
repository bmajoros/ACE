/****************************************************************
 ACE.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_ACE_H
#define INCL_ACE_H

class ACE {
public:
  ACE();
  ~ACE();
private:
  FastaWriter fastaWriter;
  SignalSensors sensors;
  String labelingFile, refProteinFile, altProteinFile;
  String substrate, altDefline, xmlFilename;
  int MAX_SPLICE_SHIFT, MIN_EXON_LEN, MIN_INTRON_LEN;
  bool allowExonSkipping, allowIntronRetention, allowCrypticSites;
  bool reverseCigar;
  String CIGAR;
  Essex::CompositeNode *root, *startCodonMsg;
  Regex warningsRegex, errorsRegex;
  int VCFwarnings, VCFerrors;
  GarbageIgnorer garbageCollector;
  GffTranscript *loadGff(const String &filename);
  void parseNoncanonicals(const String &,Set<String> &);
  String loadSeq(const String &filename);
  String loadSeq(const String &filename,String &cigar);
  void computeLabeling(GffTranscript &,Labeling &);
  void mapLabeling(Labeling &from,Labeling &to,const CigarString &);
  void mapTranscript(GffTranscript &,const CigarString &,
		     const String &outfile,const String &genome,
		     const Sequence &genomeSeq);
  bool mapExon(GffExon &,CigarAlignment &);
  void writeProtein(const String &defline,const String &protein,
		    const String &filename);
  void append(Essex::CompositeNode *,const String &tag,const String &message);
  void append(Essex::CompositeNode *,const char *tag,const char *message);
  void append(Essex::CompositeNode *,const char *tag,const String &message);
  void append(Essex::CompositeNode *,const char *tag,int);
  void appendBrokenSignals(const TranscriptSignals *,
			   Essex::CompositeNode *status);
  void writeXML();
  void processConfig(const String &filename);
  void parseConsensusList(const String &tag,ConfigFile &,Set<String> &into);
  SignalSensor *loadModel(const String &label,ConfigFile &);
};

#endif

