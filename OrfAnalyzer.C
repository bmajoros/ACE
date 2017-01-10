/****************************************************************
 OrfAnalyzer.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "OrfAnalyzer.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/CodonIterator.H"
#include "BOOM/Constants.H"
#include "SignalPrinter.H"
using namespace std;
using namespace BOOM;


OrfAnalyzer::OrfAnalyzer(SignalSensors &sensors,int MIN_ORF_LEN)
  : sensors(sensors), MIN_ORF_LEN(MIN_ORF_LEN)
{
}



GffTranscript *OrfAnalyzer::findORF(const GffTranscript &original,
				    const String &genomeStr,
				    const Sequence &genomeSeq,
				    double &startCodonScore,
				    int &genomicStartPos,
				    int &orfLength)
{
  GffTranscript *transcript=new GffTranscript(original);
  transcript->loadSequence(genomeStr);
  //String RNA=original.getFullSequence();
  String RNA=transcript->getFullSequence();
  int splicedStartPos=findStartCodon(RNA,startCodonScore);
  if(splicedStartPos<0) { delete transcript; return NULL; }
  genomicStartPos=transcript->mapToGenomicCoords(splicedStartPos);
  transcript->forgetCDS();
  transcript->splitUTRandCDS(genomeStr,genomicStartPos,sensors.stopCodons);
  orfLength=transcript->getCDSlength();
  return transcript;
}



double OrfAnalyzer::startCodonScore(const String &RNA,int consensusPos,
				   SignalSensor *&sensor)
{
  sensor=sensors.startCodonSensor;
  const int footprint=sensor->getContextWindowLength();
  const int offset=sensor->getConsensusOffset();
  const int consensusLength=sensor->getConsensusLength();
  const int L=RNA.length();
  Sequence seq(RNA,DnaAlphabet::global());
  const int last=L-footprint;
  int start=consensusPos-offset;
  int end=start+footprint;
  if(start<0 || start>last) {
    sensor=sensors.shortStartSensor;
    start=consensusPos; end=consensusPos+consensusLength;
    if(start<0 || start>last) return NEGATIVE_INFINITY;
  }
  if(!sensor->consensusOccursAt(RNA,consensusPos)) return NEGATIVE_INFINITY;
  double score=sensor->getLogP(seq,RNA,start);
  return score;
}



int OrfAnalyzer::findStartCodon(const String &transcript,
				double &startCodonScore)
{
  SignalSensor *sensor=sensors.startCodonSensor;
  const int footprint=sensor->getContextWindowLength();
  const int offset=sensor->getConsensusOffset();
  const int L=transcript.length();
  const double cutoff=sensor->getCutoff();
  Sequence seq(transcript,DnaAlphabet::global());
  const int last=L-footprint;
  for(int pos=0 ; pos<last ; ++pos) {
    if(!sensor->consensusOccursAt(transcript,pos+offset)) continue;
    startCodonScore=sensor->getLogP(seq,transcript,pos);
    if(startCodonScore>=cutoff) return pos+offset;
  }
  return -1;
}



bool OrfAnalyzer::findStartCodons(const String &transcript,
				  Vector<StartCodon> &hits)
{
  SignalSensor *sensor=sensors.startCodonSensor;
  const int footprint=sensor->getContextWindowLength();
  const int offset=sensor->getConsensusOffset();
  const int L=transcript.length();
  const double cutoff=sensor->getCutoff();
  Sequence seq(transcript,DnaAlphabet::global());
  const int last=L-footprint;
  for(int pos=0 ; pos<last ; ++pos) {
    if(!sensor->consensusOccursAt(transcript,pos+offset)) continue;
    double score=sensor->getLogP(seq,transcript,pos);
    if(score>=cutoff) hits.push_back(StartCodon(pos+offset,score));
  }
  return hits.size()>0;
}



Essex::CompositeNode *OrfAnalyzer::noncodingToCoding(
				    const GffTranscript &refTrans,
				    const String &refStr,
				    const Sequence &refSeq,
				    const GffTranscript &altTrans,
				    const String &altStr,
				    const Sequence &altSeq,
				    int &refOrfLen,
				    int &altOrfLen,
				    bool reverseStrand,
				    int altSeqLen,
				    double &refStartScore,
				    double &altStartScore,
				    String &msg)
{
  refStartScore=double(NEGATIVE_INFINITY);
  altStartScore=double(NEGATIVE_INFINITY);
  refOrfLen=altOrfLen=0;
  int refGenomicStart, altGenomicStart;
  GffTranscript *refORF=findORF(refTrans,refStr,refSeq,refStartScore,
				refGenomicStart,refOrfLen);
  GffTranscript *altORF=findORF(altTrans,altStr,altSeq,altStartScore,
				altGenomicStart,altOrfLen);
  bool change=false;
  if(!refORF && altORF && altOrfLen>=MIN_ORF_LEN)
    { change=true; msg="ref-no-start-codon"; }
  else if(refORF && altORF && refOrfLen<MIN_ORF_LEN && 
	  altOrfLen>=MIN_ORF_LEN && altOrfLen>=2*refOrfLen)
    { change=true; msg="ref-ORF-too-short"; }
  else if(refORF && altORF && altOrfLen>=MIN_ORF_LEN && altOrfLen>=refOrfLen)
    { change=true; msg="possible-misannotation"; }
  if(change) {
    altORF->computePhases();
    altORF->loadSequence(altStr);
    if(reverseStrand) altORF->reverseComplement(altSeqLen);
  }
  Essex::CompositeNode *ret=change ? altORF->toEssex() : NULL;
  delete refORF; delete altORF;
  return ret;
}



GffTranscript *
OrfAnalyzer::earlierStartCodon(const GffTranscript &refTrans,
			       const String &refStr,
			       const Sequence &refSeq,
			       const GffTranscript &altTrans,
			       const String &altStr,
			       const Sequence &altSeq,
			       const CigarAlignment &altToRef,
			       int &oldOrfLen,
			       int &newOrfLen,
			       double &oldStartCodonScore,
			       double &newStartCodonScore,
			       String &oldStartStr,
			       String &newStartStr,
			       bool reverseStrand,
			       int altSeqLen,
			       Essex::CompositeNode *&msg)
{
  /* Requirements: either the new start codon didn't exist in the reference,
     or it had a much weaker score, or it was in a different reading frame.
   */

  msg=NULL;

  // See if there is an upstream start codon
  int altGenomicStart;
  GffTranscript *altORF=findORF(altTrans,altStr,altSeq,newStartCodonScore,
				altGenomicStart,newOrfLen);
  if(!altORF) return NULL;
  oldOrfLen=altTrans.getCDSlength();
  int oldBegin, oldEnd;
  altTrans.getCDSbeginEnd(oldBegin,oldEnd);
  if(altGenomicStart>=oldBegin) { delete altORF; return NULL; }
  SignalSensor *sensor=sensors.startCodonSensor;

  // If that upstream start codon was previously intronic, it's a good bet
  // that this is a functional change
  const int refGenomicStart=altToRef.mapApproximate(altGenomicStart,DIR_RIGHT);
  const int refLocalStart=refTrans.mapToTranscriptCoords(refGenomicStart);
  String reason;
  bool change=refLocalStart<0; // -1 means unmapped due to being intronic
  if(refLocalStart<0) reason="was-intronic";

  // If this start codon existed but had a poor score in the reference,
  // it may indeed indicate a functional change
  const int offset=sensor->getConsensusOffset();
  const int windowLen=sensor->getContextWindowLength();
  if(refLocalStart>=0) { // it was exonic in the ref annotation
    GffTranscript refCopy(refTrans);
    refCopy.loadSequence(refStr);
    String refRNA=refCopy.getFullSequence();
    const int begin=refLocalStart-offset;
    if(!sensor->consensusOccursAt(refRNA,refLocalStart)) {
      change=true;
      reason="bad-consensus";
    }
    else if(begin>=0) {
      Sequence rnaSeq(refRNA,DnaAlphabet::global());
      double refScore=sensor->getLogP(rnaSeq,refRNA,begin);
      if(refScore<sensor->getCutoff()) {
	change=true;
	reason="score-below-threshold";
      }
    }
  }

  // If no change, just return
  if(!change) return NULL;

  // Compute oldStartCodonScore 
  const int altLocal=altTrans.mapToTranscriptCoords(oldBegin);
  GffTranscript altCopy(altTrans);
  altCopy.loadSequence(altStr);
  String altRNA=altCopy.getFullSequence();
  Sequence altRNAseq(altRNA,DnaAlphabet::global());
  oldStartCodonScore=sensor->getLogP(altRNAseq,altRNA,altLocal-offset);
  oldStartStr=SignalPrinter::print(*sensor,altLocal-offset,altRNA);
  const int newAltLocal=altTrans.mapToTranscriptCoords(altGenomicStart);
  newStartStr=SignalPrinter::print(*sensor,newAltLocal-offset,altRNA);

  msg=new Essex::CompositeNode("reference");
  msg->append(reason);
  altORF->computePhases();
  altORF->loadSequence(altStr);
  return altORF;
}



Essex::CompositeNode *OrfAnalyzer::lostUORFs(const GffTranscript &refTrans,
					     const String &refStr,
					     const Sequence &refSeq,
					     const GffTranscript &altTrans,
					     const String &altStr,
					     const Sequence &altSeq,
					     const CigarAlignment &refToAlt,
					     bool reverseStrand)
{
  GffTranscript refCopy(refTrans), altCopy(altTrans);
  refCopy.loadSequence(refStr); altCopy.loadSequence(altStr);
  String refRNA=refCopy.getFullSequence(), altRNA=altCopy.getFullSequence();
  Vector<StartCodon> refStarts; // spliced coordinates
  bool found=findStartCodons(refRNA,refStarts);
  if(!found) return NULL;
  Essex::CompositeNode *root=new Essex::CompositeNode("uORFs");
  Vector<GffExon*> refExons, altExons;
  refTrans.getRawExons(refExons); altTrans.getRawExons(altExons);
  GffTranscript::sort(refExons); GffTranscript::sort(altExons);
  int refStart, refEnd;
  refTrans.getCDSbeginEnd(refStart,refEnd);
  refStart= GffTranscript::genomicToSplicedCoords(refStart,refExons);
  for(Vector<StartCodon>::iterator cur=refStarts.begin(), end=refStarts.end() ;
      cur!=end ; ++cur) {
    StartCodon codon=*cur;
    if(codon.pos>=refStart) continue; // main ORF
    const int splicedRefPos=codon.pos;
    int genomicRefPos=
      GffTranscript::splicedToGenomicCoords(splicedRefPos,refExons);
    int genomicAltPos=refToAlt.mapApproximate(genomicRefPos);
    if(genomicAltPos<0) INTERNAL_ERROR;
    const int splicedAltPos=
      GffTranscript::genomicToSplicedCoords(genomicAltPos,altExons);
    double altScore=NEGATIVE_INFINITY;
    double cutoff=sensors.startCodonSensor->getCutoff();
    if(splicedAltPos>=0) {
      SignalSensor *sensor;
      altScore=startCodonScore(altRNA,splicedAltPos,sensor);
      cutoff=sensor->getCutoff();
    }
    if(isFinite(altScore) && altScore>=cutoff) {
      /*
      Essex::CompositeNode *uORFnode=new Essex::CompositeNode("retained-uORF");
      Essex::CompositeNode *refNode=new Essex::CompositeNode("ref");
      refNode->append(genomicRefPos); refNode->append(codon.score);
      uORFnode->append(refNode);
      Essex::CompositeNode *altNode=new Essex::CompositeNode("alt");
      altNode->append(genomicAltPos); altNode->append(altScore);
      uORFnode->append(altNode);
      root->append(uORFnode);
      */
    }
    else {
      int orfLen;
      {
	GffTranscript temp(refTrans);
	temp.forgetCDS();
	temp.splitUTRandCDS(refStr,genomicRefPos,sensors.stopCodons);
	orfLen=temp.getCDSlength();
      }
      if(reverseStrand) {
	genomicRefPos=refStr.length()-genomicRefPos;
	genomicAltPos=altStr.length()-genomicAltPos; }
      Essex::CompositeNode *uORFnode=new Essex::CompositeNode("lost-uORF");
      Essex::CompositeNode *refNode=new Essex::CompositeNode("ref");
      refNode->append(genomicRefPos); refNode->append(codon.score);
      refNode->append("length:"); refNode->append(orfLen);
      uORFnode->append(refNode);
      Essex::CompositeNode *altNode=new Essex::CompositeNode("alt");
      altNode->append(genomicAltPos); altNode->append(altScore);
      SignalSensor *sensor=sensors.startCodonSensor;
      int refWindow=splicedRefPos-sensor->getConsensusOffset();
      int altWindow=splicedAltPos-sensor->getConsensusOffset();
      if(altWindow<0 || 
	 altWindow+sensor->getContextWindowLength()>altRNA.length()) {
	sensor=sensors.shortStartSensor;
	refWindow=splicedRefPos; 
	altWindow=splicedAltPos; }
      String oldSignal=SignalPrinter::print(*sensor,refWindow,refRNA);
      //cout<<splicedAltPos<<" "<<sensor->getSignalType()<<" "<<altWindow<<" "<<altRNA.length()<<endl;
      String newSignal=splicedAltPos>=0 ? 
	SignalPrinter::print(*sensor,altWindow,altRNA) : "intronic";
      refNode->append(oldSignal); altNode->append(newSignal);
      altNode->append("threshold:");
      altNode->append(double(sensor->getCutoff()));
      uORFnode->append(altNode);
      root->append(uORFnode);
    }
  }
  GffTranscript::deleteExons(refExons); GffTranscript::deleteExons(altExons);
  if(root->getNumChildren()>0) return root;
  delete root;
  return NULL;
}



