/****************************************************************
 ProjectionChecker.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/Constants.H"
#include "ProjectionChecker.H"
#include "NMD.H"
using namespace std;
using namespace BOOM;



ProjectionChecker::ProjectionChecker(GffTranscript &refTrans,
				     GffTranscript &altTrans,
				     const String &refSubstrate,
				     const Sequence &refSeq,
				     const String &altSubstrate,
				     const Sequence &altSeq,
				     const Labeling &labeling,
				     const SignalSensors &sensors)
  : refTrans(refTrans), altTrans(altTrans),
    refSubstrate(refSubstrate), altSubstrate(altSubstrate),
    labeling(labeling), refSeq(refSeq), altSeq(altSeq), sensors(sensors)
{
  // ctor
}



bool ProjectionChecker::checkSpliceSiteStrengths()
{
  const int numExons=refTrans.getNumExons();
  bool ok=true;
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &refExon=refTrans.getIthExon(i);
    GffExon &altExon=altTrans.getIthExon(i);
    if(refExon.hasDonor()) 
      ok=checkDonorStrength(refExon,refSubstrate,altExon,altSubstrate) 
	&& ok;
    if(refExon.hasAcceptor()) 
      ok=checkAcceptorStrength(refExon,refSubstrate,altExon,altSubstrate)
	&& ok;
  }
  const int numUTR=refTrans.numUTR();
  for(int i=0 ; i<numUTR ; ++i) {
    GffExon &refExon=refTrans.getIthUTR(i);
    GffExon &altExon=altTrans.getIthUTR(i);
    if(refExon.hasDonor()) 
      ok=checkDonorStrength(refExon,refSubstrate,altExon,altSubstrate)
	&& ok;
    if(refExon.hasAcceptor()) 
      ok=checkAcceptorStrength(refExon,refSubstrate,altExon,altSubstrate)
	&& ok;
  }
  return ok;
}



float ProjectionChecker::scoreDonor(GffExon &exon,const String &str,
				    const Sequence &seq)
{
  if(exon.getStrand()!=FORWARD_STRAND) 
    throw "Reverse-strand features are not yet supported";
  const int donorPos=exon.getEnd();
  const int windowPos=donorPos-sensors.donorSensor->getConsensusOffset();
  if(windowPos<0 || 
     windowPos+sensors.donorSensor->getContextWindowLength()>str.length())
    return NEGATIVE_INFINITY;
  const float score=sensors.donorSensor->getLogP(seq,str,windowPos);
  return score;
}



float ProjectionChecker::scoreAcceptor(GffExon &exon,const String &str,
				       const Sequence &seq)
{
  if(exon.getStrand()!=FORWARD_STRAND) 
    throw "Reverse-strand features are not yet supported";
  const int acceptorPos=exon.getBegin()-2;
  const int windowPos=acceptorPos-sensors.acceptorSensor->getConsensusOffset();
  if(windowPos<0 || 
     windowPos+sensors.acceptorSensor->getContextWindowLength()>str.length())
    return NEGATIVE_INFINITY;
  const float score=sensors.acceptorSensor->getLogP(seq,str,windowPos);
  return score;
}



bool ProjectionChecker::checkDonorStrength(GffExon &refExon,
					   const String &refStr,
					   GffExon &altExon,
					   const String &altStr)
{
  const float refScore=scoreDonor(refExon,refStr,refSeq);
  if(refScore<sensors.donorSensor->getCutoff()) return true;
  const float altScore=scoreDonor(altExon,altStr,altSeq);
  if(altScore<sensors.donorSensor->getCutoff()) {
    cout<<"Donor splice site has been weakened: "<<refScore
	<<" ("<<refExon.getEnd()<<" in ref) vs. "<<altScore
	<<" ("<<altExon.getEnd()<<" in alt)"<<endl;
    return false;
  }
  return true;
}
  


bool ProjectionChecker::checkAcceptorStrength(GffExon &refExon,
					      const String &refStr,
					      GffExon &altExon,
					      const String &altStr)
{
  const float refScore=scoreAcceptor(refExon,refStr,refSeq);
  if(refScore<sensors.acceptorSensor->getCutoff()) return true;
  const float altScore=scoreAcceptor(altExon,altStr,altSeq);
  if(altScore<sensors.acceptorSensor->getCutoff()) {
    cout<<"Acceptor splice site has been weakened: "<<refScore
	<<" ("<<refExon.getBegin()-2<<" in ref) vs. "<<altScore
	<<" ("<<altExon.getBegin()-2<<" in alt)"<<endl;
    return false;
  }
  return true;
}



TranscriptSignals *ProjectionChecker::findBrokenSpliceSites()
{
  Vector<GffExon*> refExons, altExons;
  refTrans.getRawExons(refExons);
  altTrans.getRawExons(altExons);
  if(altTrans.getStrand()!=FORWARD_STRAND) INTERNAL_ERROR;
  TranscriptSignals *signals=NULL;
  if(refExons.size()==altExons.size()) {
    const int numExons=refExons.size();
    signals=new TranscriptSignals;
    signals->setID(refTrans.getTranscriptId());
    signals->setSubstrate(altTrans.getSubstrate());
    signals->setSource("ICE");
    signals->setStrand(refTrans.getStrand());
    signals->setGeneID(refTrans.getGeneId());
    for(int i=0 ; i<numExons ; ++i) {
      GffExon &refExon=*refExons[i];
      GffExon &altExon=*altExons[i];
      const int altBegin=altExon.getBegin(), altEnd=altExon.getEnd();
      if(i==0) signals->addSignal(TSS,altBegin,0.0);
      else { // i>0
	bool weakened; String consensus, window;
	float refScore, altScore, cutoff;
	bool broken=!checkAcceptor(refExon,altExon,weakened,consensus,window,
				   refScore,altScore,cutoff);
	TranscriptSignal &signal=signals->addSignal(AG,altBegin-2,altScore);
	signal.broken=broken; signal.weakened=weakened;
	signal.seq=window ;signal.refScore=refScore; signal.cutoff=cutoff;
      }
      if(i<numExons-1) {
	bool weakened; String consensus, window;
	float refScore, altScore, cutoff;
	bool broken=!checkDonor(refExon,altExon,weakened,consensus,window,
				refScore,altScore,cutoff);
	TranscriptSignal &signal=signals->addSignal(GT,altEnd,altScore);
	signal.broken=broken; signal.weakened=weakened;
	signal.seq=window; signal.refScore=refScore; signal.cutoff=cutoff;
      }
      else signals->addSignal(TES,altEnd,0.0);
    }
    bool isCoding=altTrans.numExons()>0; // actually # coding segments
    if(isCoding) {
      int start=altTrans.getIthExon(0).getBegin(); // assumes forward strand
      signals->setStartCodon(start);
    }
  }
  GffTranscript::deleteExons(refExons);
  GffTranscript::deleteExons(altExons);
  return signals;
}



/*
TranscriptSignals *ProjectionChecker::simulateBrokenSpliceSites()
{
  Vector<GffExon*> refExons, altExons;
  refTrans.getRawExons(refExons);
  altTrans.getRawExons(altExons);
  if(altTrans.getStrand()!=FORWARD_STRAND) INTERNAL_ERROR;
  TranscriptSignals *signals=NULL;
  if(refExons.size()==altExons.size()) {
    const int numExons=refExons.size();
    signals=new TranscriptSignals;
    signals->setID(refTrans.getTranscriptId());
    signals->setSubstrate(altTrans.getSubstrate());
    signals->setSource("SIMULATION");
    throw "simulation";
    signals->setStrand(refTrans.getStrand());
    signals->setGeneID(refTrans.getGeneId());
    for(int i=0 ; i<numExons ; ++i) {
      GffExon &refExon=*refExons[i];
      GffExon &altExon=*altExons[i];
      const int altBegin=altExon.getBegin(), altEnd=altExon.getEnd();
      if(i==0) signals->addSignal(TSS,altBegin,0.0);
      else { // i>0
	bool weakened; String consensus, window;
	float refScore, altScore, cutoff;
	bool broken=!checkAcceptor(refExon,altExon,weakened,consensus,window,
				   refScore,altScore,cutoff);
	TranscriptSignal &signal=signals->addSignal(AG,altBegin-2,altScore);
	signal.broken=broken; signal.weakened=weakened;
	signal.seq=window ;signal.refScore=refScore; signal.cutoff=cutoff;
      }
      if(i<numExons-1) {
	bool weakened; String consensus, window;
	float refScore, altScore, cutoff;
	bool broken=!checkDonor(refExon,altExon,weakened,consensus,window,
				refScore,altScore,cutoff);
	TranscriptSignal &signal=signals->addSignal(GT,altEnd,altScore);
	signal.broken=broken; signal.weakened=weakened;
	signal.seq=window; signal.refScore=refScore; signal.cutoff=cutoff;
      }
      else signals->addSignal(TES,altEnd,0.0);
    }
    bool isCoding=altTrans.numExons()>0; // actually # coding segments
    if(isCoding) {
      int start=altTrans.getIthExon(0).getBegin(); // assumes forward strand
      signals->setStartCodon(start);
    }
  }
  GffTranscript::deleteExons(refExons);
  GffTranscript::deleteExons(altExons);
  if(signals) signals->simulateBroken();
  return signals;
}
*/



bool ProjectionChecker::checkSpliceSites(bool quiet)
{
  const int numExons=refTrans.getNumExons();
  if(altTrans.getNumExons()!=numExons) {
    if(!quiet) 
      cout<<"Projected transcript has different number of exons"<<endl;
    return false;
  }
  bool ok=true;
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &refExon=refTrans.getIthExon(i);
    GffExon &altExon=altTrans.getIthExon(i);
    bool weakened; String consensus, window;
    float refScore, altScore, cutoff;
    if(refExon.hasDonor()) 
      ok=checkDonor(refExon,altExon,weakened,consensus,window,
		    refScore,altScore,cutoff) && ok;
    if(refExon.hasAcceptor()) 
      ok=checkAcceptor(refExon,altExon,weakened,consensus,window,
		       refScore,altScore,cutoff) && ok;
  }
  const int numUTR=refTrans.numUTR();
  for(int i=0 ; i<numUTR ; ++i) {
    GffExon &refExon=refTrans.getIthUTR(i);
    GffExon &altExon=altTrans.getIthUTR(i);
    bool weakened; String consensus, window;
    float refScore, altScore, cutoff;
    if(refExon.hasDonor()) 
      ok=checkDonor(refExon,altExon,weakened,consensus,window,
		    refScore,altScore,cutoff) && ok;
    if(refExon.hasAcceptor()) 
      ok=checkAcceptor(refExon,altExon,weakened,consensus,window,
		       refScore,altScore,cutoff) && ok;
  }

  return ok;
}



String ProjectionChecker::getParsedWindow(SignalSensor &sensor,
					  int windowPos,const String &genome)
{
  String window;
  const int offset=sensor.getConsensusOffset();
  const int windowLen=sensor.getContextWindowLength();
  window=
    genome.substring(windowPos,offset).tolower()
    + "_"
    + genome.substring(windowPos+offset,2) 
    + "_"
    + genome.substring(windowPos+offset+2,windowLen-offset-2).tolower();
  return window;
}



String ProjectionChecker::spliceSiteChangeString(SignalSensor *sensor,
						 int refBegin,
						 const String &refSubstrate,
						 float refScore,
						 int altBegin,
						 const String &altSubstrate,
						 float altScore)
{
  return
    getParsedWindow(*sensor,refBegin,refSubstrate)
    + " " + refScore + " "
    + getParsedWindow(*sensor,altBegin,altSubstrate)
    + " " + altScore;
}



bool ProjectionChecker::checkDonor(GffExon &refExon,GffExon &altExon,
				   bool &weakened,String &altDonor,
				   String &altWindow,float &refScore,
				   float &altScore,float &cutoff)
{
  // NOTE: this function assumes forward-strand features only ###

  weakened=false;
  refScore=scoreDonor(refExon,refSubstrate,refSeq);
  altScore=scoreDonor(altExon,altSubstrate,altSeq);
  cutoff=sensors.donorSensor->getCutoff();
  int refPos, pos;
  const String refDonor=getDonor(refExon,refSubstrate,refPos);
  altDonor=getDonor(altExon,altSubstrate,pos);
  if(altDonor!=refDonor && !sensors.donorConsensuses.isMember(altDonor)) {
    const int offset=sensors.donorSensor->getConsensusOffset();
    const int altBegin=altExon.getEnd()-offset;
    const int refBegin=refExon.getEnd()-offset;
    altWindow=spliceSiteChangeString(sensors.donorSensor,refBegin,refSubstrate,
				     refScore,altBegin,altSubstrate,altScore);
    return false;
  }

  // The reference must score above threshold, the alt must score below
  // threshold, and the difference must be at least a factor of 2:
  if(refScore>=cutoff && altScore<cutoff && refScore-altScore>=log(2)) { 
    const int offset=sensors.donorSensor->getConsensusOffset();
    const int altBegin=altExon.getEnd()-offset;
    const int refBegin=refExon.getEnd()-offset;
    altWindow=spliceSiteChangeString(sensors.donorSensor,refBegin,refSubstrate,
				     refScore,altBegin,altSubstrate,altScore);
    weakened=true;
    //return false;
  }
  return true;
}  



bool ProjectionChecker::checkAcceptor(GffExon &refExon,GffExon &altExon,
				      bool &weakened,String &altAcceptor,
				      String &altWindow,float &refScore,
				      float &altScore,float &cutoff)
{
  // NOTE: this function assumes forward-strand features only ###

  weakened=false;
  refScore=scoreAcceptor(refExon,refSubstrate,refSeq);
  altScore=scoreAcceptor(altExon,altSubstrate,altSeq);
  cutoff=sensors.acceptorSensor->getCutoff();
  int refPos, pos;
  const String refAcceptor=getAcceptor(refExon,refSubstrate,refPos);
  altAcceptor=getAcceptor(altExon,altSubstrate,pos);
  if(altAcceptor!=refAcceptor && 
     !sensors.acceptorConsensuses.isMember(altAcceptor)) {
    const int offset=sensors.acceptorSensor->getConsensusOffset();
    const int altBegin=altExon.getBegin()-2-offset;
    const int refBegin=refExon.getBegin()-2-offset;
    altWindow=spliceSiteChangeString(sensors.acceptorSensor,refBegin,
				     refSubstrate,refScore,altBegin,
				     altSubstrate,altScore);
    return false;
    }
  if(refScore>=cutoff && altScore<cutoff && refScore-altScore>=log(2)) { 
    const int offset=sensors.acceptorSensor->getConsensusOffset();
    const int altBegin=altExon.getBegin()-2-offset;
    const int refBegin=refExon.getBegin()-2-offset;
    altWindow=spliceSiteChangeString(sensors.acceptorSensor,refBegin,
				     refSubstrate,refScore,altBegin,
				     altSubstrate,altScore);
    weakened=true; 
    //return false; 
  }
  return true;
}



String ProjectionChecker::getDonor(GffExon &exon,const String &substrate,
				   int &pos)
{
  if(exon.getStrand()=='+') {
    const int end=exon.getEnd();
    if(end>substrate.length()-2) return "";
    return substrate.substring(pos=end,2);
  }
  else {
    const int begin=exon.getBegin();
    if(begin<2 || begin>substrate.length()) return "";
    return ProteinTrans::reverseComplement(substrate.substring(pos=begin-2,2));
  }
}



String ProjectionChecker::getAcceptor(GffExon &exon,const String &substrate,int &pos)
{
  if(exon.getStrand()=='+') {
    const int begin=exon.getBegin();
    if(begin<2 || begin>substrate.length()) return "";
    return substrate.substring(pos=begin-2,2);
  }
  else {
    const int end=exon.getEnd();
    if(end+2>substrate.length()) return "";
    return ProteinTrans::reverseComplement(substrate.substring(pos=end,2));
  }
}



bool ProjectionChecker::checkFrameshifts(const Labeling &labeling,
					 const GffTranscript &transcript,
					 const String &substrate,
					 Essex::CompositeNode *status)
{
  if(labeling.length()!=substrate.length()) {
    status->append("internal-error",
		   "labeling and alt substrate have different lengths");
    return false;
  }
  const int numExons=transcript.numExons();
  int phase=0, phaseMatches=0, phaseMismatches=0;
  for(int i=0 ; i<numExons ; ++i) {
    const GffExon &exon=transcript.getIthExon(i);
    const int begin=exon.getBegin(), end=exon.getEnd();
    for(int pos=begin ; pos<end ; ++pos) {
      const GeneModelLabel label=labeling[pos];
      if(isExon(label))
	if(phase==getExonPhase(label)) ++phaseMatches;
	else ++phaseMismatches;
      phase=(phase+1)%3;
    }
  }
  if(phaseMismatches>0) {
    const int total=phaseMismatches+phaseMatches;
    float percentMismatch=int(1000*phaseMismatches/float(total)+5/9.0)/10.0;
    Essex::CompositeNode *node=new Essex::CompositeNode("frameshift");
    node->append("nt-with-phase-mismatch",phaseMismatches);
    node->append("percent-phase-mismatch",String(percentMismatch)+"%");
    status->append(node);
    return false;
  }
  return true;
}



bool ProjectionChecker::checkFrameshifts(const Labeling &refLab,
					 const Labeling &altLab,
					 Essex::CompositeNode *status)
{
  if(refLab.length()!=altLab.length()) {
    status->append("internal-error",
		   "ref and alt labelings have different lengths");
    return false;
  }
  const int L=refLab.length();
  int phaseMismatches=0, phaseMatches=0;
  for(int i=0 ; i<L ; ++i) {
    const GeneModelLabel ref=refLab[i];
    if(isExon(ref)) {
      GeneModelLabel alt=altLab[i];
      if(isExon(alt))
	if(getExonPhase(alt)==getExonPhase(ref)) ++phaseMatches;
	else ++phaseMismatches;
    }
  }
  if(phaseMismatches>0) {
    const int total=phaseMismatches+phaseMatches;
    float percentMismatch=int(1000*phaseMismatches/float(total)+5/9.0)/10.0;
    Essex::CompositeNode *node=new Essex::CompositeNode("frameshift");
    node->append("nt-with-phase-mismatch",phaseMismatches);
    node->append("percent-phase-mismatch",String(percentMismatch)+"%");
    status->append(node);
    return false;
  }
  return true;
}



void ProjectionChecker::translate(GffTranscript &refTrans,
				  GffTranscript &altTrans,
				  String &refProtein,
				  String &altProtein)
{
  refTrans.loadSequence(refSubstrate); 
  String refDNA=refTrans.getSequence();
  refProtein=ProteinTrans::translate(refDNA);
  if(refProtein.lastChar()!='*') {
    refTrans.extendFinalExonBy3(); altTrans.extendFinalExonBy3();
    refTrans.loadSequence(refSubstrate); 
    refDNA=refTrans.getSequence();
    refProtein=ProteinTrans::translate(refDNA);
  }
  altTrans.loadSequence(altSubstrate);
  const String altDNA=altTrans.getSequence();
  altProtein=ProteinTrans::translate(altDNA);
}



bool ProjectionChecker::hasStartCodon(const String &protein)
{
  return protein.length()>0 && protein[0]=='M';
}



bool ProjectionChecker::hasStopCodon(const String &protein)
{
  return protein.lastChar()=='*';
}



bool ProjectionChecker::hasPTC(const String &protein,int &PTCpos)
{
  PTCpos=protein.findFirst('*');
  return PTCpos>=0 && PTCpos<protein.length()-1;
}



bool ProjectionChecker::geneIsWellFormed(GffTranscript &transcript,
					 const String &substrate,
					 bool &noStart,bool &noStop,
					 bool &PTC,bool &badSpliceSite,
					 Essex::CompositeNode *status,
					 const SignalSensors &sensors,
					 const int nmdDistParm)
{
  noStart=noStop=PTC=badSpliceSite=false;
  transcript.loadSequence(substrate);
  NMD nmd(nmdDistParm);
  int ejcDistance;
  NMD_TYPE nmdType=nmd.predict(transcript,substrate,ejcDistance);
  if(nmdType==NMD_NMD) {
    Essex::CompositeNode *msg=new Essex::CompositeNode("NMD");
    msg->append("EJC-distance",ejcDistance);
    status->append(msg);
    PTC=true;
    return false;
  }
  if(transcript.isCoding()) {
    String protein=ProteinTrans::translate(transcript.getSequence());
    if(protein.length()<1) {
      noStart=noStop=true;
      return false; }
    if(protein[0]!='M') {
      if(transcript.getStrand()!=FORWARD_STRAND) 
	throw "reverse strand not implemented";
      int begin=transcript.getIthExon(0).getBegin();
      if(begin>3 && substrate.substr(begin-3,3)=="ATG")
	transcript.getIthExon(0).setBegin(begin-3);
      else {
	status->append(new Essex::CompositeNode("bad-start"));
	noStart=true; }
    }
    const int firstStop=protein.findFirst('*');
    const int proteinLen=protein.length();
    if(firstStop>=0 && firstStop<proteinLen-1) {
      Essex::CompositeNode *msg=new Essex::CompositeNode("premature-stop");
      msg->append("AA-pos",firstStop);
      msg->append("protein-length",proteinLen);
      status->append(msg);
      PTC=true; }
    else if(protein.lastChar()!='*') {
      transcript.extendFinalExonBy3();
      transcript.loadSequence(substrate);
      String protein=ProteinTrans::translate(transcript.getSequence());
      if(protein.lastChar()!='*') {
	status->append(new Essex::CompositeNode("no-stop-codon"));
	noStop=true; }
    }
  }
  const int numExons=transcript.getNumExons();
  String msg;
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript.getIthExon(i);
    if(exon.hasDonor()) 
      if(!checkDonor(exon,substrate,status,sensors)) badSpliceSite=true;
    if(exon.hasAcceptor()) 
      if(!checkAcceptor(exon,substrate,status,sensors)) badSpliceSite=true;
  }
  const int numUTR=transcript.numUTR();
  for(int i=0 ; i<numUTR ; ++i) {
    GffExon &exon=transcript.getIthUTR(i);
    if(exon.hasDonor()) 
      if(!checkDonor(exon,substrate,status,sensors)) badSpliceSite=true;
    if(exon.hasAcceptor()) 
      if(!checkAcceptor(exon,substrate,status,sensors)) badSpliceSite=true;
  }
  return !(noStart || noStop || PTC || badSpliceSite);
}



bool ProjectionChecker::checkDonor(GffExon &refExon,const String &refSubstrate,
				   Essex::CompositeNode *status,
				   const SignalSensors &sensors)
{
  int pos;
  const String refDonor=getDonor(refExon,refSubstrate,pos);
  if(sensors.donorConsensuses.isMember(refDonor)) return true;
  Essex::CompositeNode *msg=new Essex::CompositeNode("bad-donor");
  msg->append(refDonor);
  msg->append(pos);
  status->append(msg);
  return false;
}



bool ProjectionChecker::checkAcceptor(GffExon &refExon,
				      const String &refSubstrate,
				      Essex::CompositeNode *status,
				      const SignalSensors &sensors)
{
  int pos;
  const String refAcceptor=getAcceptor(refExon,refSubstrate,pos);
  if(sensors.acceptorConsensuses.isMember(refAcceptor)) return true;
  Essex::CompositeNode *msg=new Essex::CompositeNode("bad-acceptor");
  msg->append(refAcceptor);
  msg->append(pos);
  status->append(msg);
  return false;
}


