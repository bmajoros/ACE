/****************************************************************
 OrfAnalyzer.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_OrfAnalyzer_H
#define INCL_OrfAnalyzer_H
#include <iostream>
#include "BOOM/GffTranscript.H"
#include "BOOM/String.H"
#include "BOOM/Sequence.H"
#include "BOOM/CigarAlignment.H"
#include "SignalSensors.H"
using namespace std;
using namespace BOOM;

struct StartCodon {
  StartCodon(int pos,float score) : pos(pos), score(score) {}
  int pos;
  float score;
};

/*  This class works only on forward-strand genes.
 */
class OrfAnalyzer {
public:
  OrfAnalyzer(SignalSensors &,int MIN_ORF_LEN);
  GffTranscript *findORF(const GffTranscript &,const String &genomeStr,
			 const Sequence &genomeSeq,float &startCodonScore,
			 int &genomicStartCodonPos,int &orfLength);
  Essex::CompositeNode *noncodingToCoding(const GffTranscript &refTrans,
					  const String &refStr,
					  const Sequence &refSeq,
					  const GffTranscript &altTrans,
					  const String &altStr,
					  const Sequence &altSeq,
					  int &refOrfLen,
					  int &altOrfLen,
					  bool reverseStrand,
					  int altSeqLen,
					  float &refStartScore,
					  float &altStartScore,
					  String &msg);
  GffTranscript *earlierStartCodon(const GffTranscript &refTrans,
				   const String &refStr,
				   const Sequence &refSeq,
				   const GffTranscript &altTrans,
				   const String &altStr,
				   const Sequence &altSeq,
				   const CigarAlignment &altToRef,
				   int &refOrfLen,
				   int &altOrfLen,
				   float &oldStartCodonScore,
				   float &newStartCodonScore,
				   String &oldStartStr,
				   String &newStartStr,
				   bool reverseStrand,
				   int altSeqLen,
				   Essex::CompositeNode *&msg);
  Essex::CompositeNode *lostUORFs(const GffTranscript &refTrans,
				  const String &refStr,
				  const Sequence &refSeq,
				  const GffTranscript &altTrans,
				  const String &altStr,
				  const Sequence &altSeq,
				  const CigarAlignment &refToAlt,
				  bool reverseStrand);
  bool findStartCodons(const String &transcript,Vector<StartCodon> &);
  float startCodonScore(const String &RNA,int consensusPos,SignalSensor *&);
private:
  SignalSensors &sensors;
  const int MIN_ORF_LEN;
  int findStartCodon(const String &transcript,float &startCodonScore);
};

#endif

