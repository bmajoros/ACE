/****************************************************************
 ACE.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "ACE.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 Globals
 ****************************************************************/
bool VERBOSE=false;
static const char *PROGRAM_NAME="ACE";
static const char *VERSION="1.0";
Alphabet alphabet;




/****************************************************************
 ACE constructor
 ****************************************************************/
ACE::ACE()
  : warningsRegex("/warnings=(\\d+)"), errorsRegex("/errors=(\\d+)"), 
    VCFwarnings(0), VCFerrors(0), startCodonMsg(NULL), substMatrix(NULL),
    variantRegex("(\\S+):(\\S+):(\\d+):(\\d+):([^:]*):([^:]*)"),
    coordRegex("/coord=(\\S+)"), orfAnalyzer(NULL),
    alignment(NULL), revAlignment(NULL), status(NULL)
{
  alphabet=DnaAlphabet::global();
}



/****************************************************************
 ACE destructor
 ****************************************************************/
ACE::~ACE()
{
  delete substMatrix;
  delete orfAnalyzer;
  delete alignment;
  delete revAlignment;
}



/****************************************************************
 ACE::main()
 ****************************************************************/
int ACE::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"ce:l:qx:v");
  parseCommandLine(cmd);

  // Read some data from files
  if(VERBOSE) cerr<<"loading inputs"<<endl;
  loadInputs(configFile,refGffFile,refFasta,altFasta);

  // Check that the reference gene is well-formed
  if(VERBOSE) cerr<<"checking reference gene"<<endl;
  status=new Essex::CompositeNode("status");
  bool referenceIsOK=checkRefGene();

  // Set up to generate structured output in Essex/XML
  if(VERBOSE) cerr<<"preparing output"<<endl;
  ofstream osACE(outACE.c_str());
  initEssex(osACE,cmd);
  
  // Compute the reference labeling
  if(VERBOSE) cerr<<"computing reference labeling"<<endl;
  Labeling refLab(refSeqLen);
  computeLabeling(*refTrans,refLab);

  // Make CIGAR alignment
  if(VERBOSE) cerr<<"building alignment"<<endl;
  buildAlignment();

  // Project the reference GFF over to an alternate GFF
  if(VERBOSE) cerr<<"mapping transcript"<<endl;
  mapTranscript(outGff);

  // Generate labeling file
  if(VERBOSE) cerr<<"projecting the labeling"<<endl;
  Labeling projectedLab(altSeqLen);
  mapLabeling(refLab,projectedLab,labelingFile);

  // Check the projection to see if the gene might be broken
  if(VERBOSE) cerr<<"checking projection"<<endl;
  bool mapped=false;
  if(referenceIsOK) 
    checkProjection(outGff,mapped,projectedLab,osACE);

  // Flush output
  if(VERBOSE) cerr<<"cleaning up"<<endl;
  flushOutput(osACE,mapped);

  cout<<"ACE terminated successfully"<<endl;
  return 0;
}



/****************************************************************
 ACE::flushOutput()
 ****************************************************************/
void ACE::flushOutput(ostream &osACE,const bool &mapped)
{
  if(mapped && status && status->getNumChildren()<2 && quiet) return;
  if(!xmlFilename.empty()) writeXML();
  osACE<<*root<<endl;
  osACE<<"#===========================================================\n";
}



/****************************************************************
 ACE::parseCommandLine()
 ****************************************************************/
void ACE::parseCommandLine(const CommandLine &cmd)
{
  if(cmd.option('v')) { cout<<"version "<<VERSION<<endl; exit(0); }
  if(cmd.numArgs()!=6)
    throw String("\n\
ace <ace.config> <ref.gff> <ref.fasta> <alt.fasta> <out.gff> <out.essex>\n\
     -c = sequence has been reversed, but cigar string has not\n\
     -e N = abort if vcf errors >N\n\
     -l <file> = emit a per-nucleotide labeling for the alt sequence\n\
     -q = quiet (only report transcripts with mapping issues)\n\
     -x <file> = also emit xml\n\
  alt.fasta must have a cigar string: >ID ... /cigar=1045M3I10M7D4023M ...\n\
\n");
  configFile=cmd.arg(0);
  refGffFile=cmd.arg(1);
  refFasta=cmd.arg(2);
  altFasta=cmd.arg(3);
  outGff=cmd.arg(4);
  outACE=cmd.arg(5);
  commandLineOpts(cmd);
}



/****************************************************************
 ACE::commandLineOpts()
 ****************************************************************/
void ACE::commandLineOpts(const CommandLine &cmd)
{
  if(cmd.option('l')) labelingFile=cmd.optParm('l');
  if(cmd.option('x')) xmlFilename=cmd.optParm('x');
  reverseStrand=cmd.option('c');
  quiet=cmd.option('q');
}



/****************************************************************
 ACE::buildAlignment()
 ****************************************************************/
void ACE::buildAlignment()
{
  cigar=CigarString(CIGAR);
  if(reverseStrand) cigar.reverse();
  alignment=cigar.getAlignment();
  revAlignment=alignment->invert(altSeqLen);
}



/****************************************************************
 ACE::loadInputs()
 ****************************************************************/
void ACE::loadInputs(const String &configFile,const String &refGffFile,
		     const String &refFasta,const String &altFasta)
{
  processConfig(configFile);
  refSeqStr=loadSeq(refFasta); altSeqStr=loadSeq(altFasta,CIGAR);
  refSeq.copyFrom(refSeqStr,alphabet); altSeq.copyFrom(altSeqStr,alphabet);
  refSeqLen=refSeqStr.length(), altSeqLen=altSeqStr.length();
  refTrans=loadGff(refGffFile);
  refTrans->loadSequence(refSeqStr);
  refProtein=refTrans->getProtein();
}



/****************************************************************
 ACE::checkRefGene()
 ****************************************************************/
bool ACE::checkRefGene()
{
  bool noStart, noStop, PTC, badSpliceSite;
  String msg;
  if(!ProjectionChecker::geneIsWellFormed(*refTrans,refSeqStr,
					  noStart,noStop,PTC,badSpliceSite,
					  status,sensors,NMD_DISTANCE_PARM)) {
    if(quiet) return false;
    status->prepend("bad-annotation");
    return false;
  }
  return true;
}



/****************************************************************
 ACE::checkProjection()
 ****************************************************************/
void ACE::checkProjection(const String &outGff,
			  bool &mapped,const Labeling &projectedLab,
			  ostream &osACE)
{
  GffTranscript *altTrans=loadGff(outGff);
  altTrans->loadSequence(altSeqStr);
  altTrans->computePhases();
  Essex::CompositeNode *altTransEssex=
    altTrans->toEssex(reverseStrand,altSeqLen);
  altTransEssex->getTag()="mapped-transcript";
  VariantClassifier classifier(variants,VariantClassifier::ALT,*altTrans,
			       altSeqLen,reverseStrand);
  altTransEssex->append(classifier.makeVariantsNode());
  root->append(altTransEssex);
  ProjectionChecker checker(*refTrans,*altTrans,refSeqStr,refSeq,
			    altSeqStr,altSeq,projectedLab,sensors);
  String altProtein;
  if(refTrans->isCoding())
    checker.translate(*refTrans,*altTrans,refProtein,altProtein);
  TranscriptSignals *signals=checker.findBrokenSpliceSites();
  if(!signals) {
    status->prepend("unequal-numbers-of-exons");
    delete altTrans;
    return;
  }

  // Enumerate alternative structures
  if(signals->anyBroken())
    { enumerateAlts(altTransEssex,signals,altTrans,osACE,projectedLab);
      delete signals; return; }
  if(signals->anyWeakened()) appendBrokenSignals(signals);
  delete signals; // ### added 7/11/2016

  // Otherwise, projection was successful
  status->prepend("mapped");
  mapped=true;

  // Translate to proteins
  if(refTrans->isCoding())
    handleCoding(altTrans,checker,projectedLab);
  else  // ref gene is noncoding
    handleNoncoding(altTrans);

  delete altTrans;
}



/****************************************************************
 ACE::initEssex()
 ****************************************************************/
void ACE::initEssex(ostream &osACE,
		    const CommandLine &cmd)
{
  String transcriptID=refTrans->getTranscriptId();
  String geneID=refTrans->getGeneId();
  root=new Essex::CompositeNode("report");
  append(root,"substrate",substrate);
  if(!globalCoord.empty()) append(root,"global-coords",globalCoord);
  append(root,"transcript-ID",transcriptID);
  append(root,"gene-ID",geneID);
  append(root,"vcf-warnings",VCFwarnings);
  append(root,"vcf-errors",VCFerrors);
  append(root,"alignment",CIGAR);
  append(root,"ref-length",refSeqStr.length());
  append(root,"alt-length",altSeqStr.length());
  //append(root,"defline",altDefline);
  //Essex::CompositeNode *essexVariants=makeEssexVariants();
  //root->append(essexVariants);
  refTrans->computePhases();
  Essex::CompositeNode *refTransEssex=
    refTrans->toEssex(reverseStrand,refSeqLen);
  refTransEssex->getTag()="reference-transcript";
  VariantClassifier classifier(variants,VariantClassifier::REF,*refTrans,
			       refSeqLen,reverseStrand);
  refTransEssex->append(classifier.makeVariantsNode());
  root->append(refTransEssex);
  root->appendChild(status);
  if(cmd.option('e') && VCFerrors>cmd.optParm('e').asInt()) {
    status->append("too-many-vcf-errors");
  }
}



Essex::Node *ACE::signalNode(const String &tag,const String &signal,
			     double score,double cutoff)
{
  Essex::CompositeNode *node=new Essex::CompositeNode(tag);
  node->append(signal);
  node->append(float(score));
  node->append("cutoff:"); node->append(float(cutoff));
  return node;
}



/****************************************************************
 ACE::analyzeEarlierStart()
 ****************************************************************/
void ACE::analyzeEarlierStart(GffTranscript *altTrans,
			      ProjectionChecker &checker,
			      const Labeling &projectedLab,
			      Essex::CompositeNode *status)
{
  // Translate
  String refProtein, altProtein;
  checker.translate(*refTrans,*altTrans,refProtein,altProtein);

  int ejcDistance;
  NMD_TYPE nmdType=nmd.predict(*altTrans,altSeqStr,ejcDistance);
  switch(nmdType) {
  case NMD_NONE: break;
  case NMD_NMD: {
    Essex::CompositeNode *stopNode=
      new Essex::CompositeNode("premature-stop");
    stopNode->append("hypothetical-NMD");
    stopNode->append("EJC-distance",ejcDistance);
    status->append(stopNode);
  } break;
  //case NMD_TRUNCATION:  { // ### this is handled below
  case NMD_NO_STOP: 
    if(refTrans->hasUTR3()) status->append("nonstop-decay");
    break;
  case NMD_NO_START: INTERNAL_ERROR;
  }

  // Check for frameshifts and amino acid differences
  const bool proteinsDiffer=refProtein!=altProtein;
  if(proteinsDiffer) {
    int matches;
    alignProteins(refProtein,altProtein,matches);
    Essex::CompositeNode *fate=new Essex::CompositeNode("protein-differs");
    percentMatch(matches,refProtein.length(),altProtein.length(),fate);
    status->append(fate);
    Labeling altLab(altSeqLen);
    computeLabeling(*altTrans,altLab);
    checker.checkFrameshifts(projectedLab,altLab,status);
  }

  // Check for premature stop codon
  if(nmdType==NMD_NONE && proteinsDiffer && refProtein.findFirst('*')>0 &&
     altProtein.findFirst('*')>0) {
    int refStop=refTrans->stopCodonGlobalCoord();
    int altStop=altTrans->stopCodonGlobalCoord();
    int altStopOnRef=(*revAlignment).mapApproximate(altStop,DIR_LEFT);
    if(altStopOnRef>=0 && refStop>=0 && altStopOnRef<refStop) {
      Essex::CompositeNode *fate=new Essex::CompositeNode("premature-stop");
      Essex::CompositeNode *truncation=
	new Essex::CompositeNode("protein-truncation");
      const int diff=getTruncationLength(*refTrans,altStopOnRef,refStop);
      truncation->append(String("")+diff+"aa");
      fate->append(truncation);
      status->append(fate);
    }
  }
}



/****************************************************************
 ACE::handleCoding()
 ****************************************************************/
void ACE::handleCoding(GffTranscript *altTrans,
		       ProjectionChecker &checker,
		       const Labeling &projectedLab)
{
  int oldOrfLen, newOrfLen; double oldStartScore, newStartScore;
  String oldStartStr, newStartStr;
  Essex::CompositeNode *reason=NULL;
  GffTranscript *newTranscript=
    orfAnalyzer->earlierStartCodon(*refTrans,refSeqStr,refSeq,
				   *altTrans,altSeqStr,altSeq,
				   *revAlignment,oldOrfLen,newOrfLen,
				   oldStartScore,newStartScore,
				   oldStartStr,newStartStr,
				   reverseStrand,altSeqLen,reason);

  if(newTranscript) {
    Essex::CompositeNode *upstreamStart;
    if(reverseStrand) {
      GffTranscript temp(*newTranscript);
      temp.reverseComplement(altSeqLen);
      upstreamStart=temp.toEssex();
    } else upstreamStart=newTranscript->toEssex();
    Essex::CompositeNode *changeNode=
      new Essex::CompositeNode("new-upstream-start-codon");
    const double cutoff=sensors.startCodonSensor->getCutoff();
    changeNode->append(signalNode("new-start",newStartStr,newStartScore,
				  cutoff));
    changeNode->append(signalNode("old-start",oldStartStr,oldStartScore,
				  cutoff));
    changeNode->append(reason);
    Essex::CompositeNode *lengthNode=new Essex::CompositeNode("ORF-length");
    lengthNode->append(oldOrfLen);
    lengthNode->append("=>");
    lengthNode->append(newOrfLen);
    changeNode->append(lengthNode);
    changeNode->append(upstreamStart);
    status->append(changeNode);
    analyzeEarlierStart(newTranscript,checker,projectedLab,changeNode);
    delete newTranscript;
  }

  // Translate
  String refProtein, altProtein;
  checker.translate(*refTrans,*altTrans,refProtein,altProtein);

  int ejcDistance;
  NMD_TYPE nmdType=nmd.predict(*altTrans,altSeqStr,ejcDistance);
  switch(nmdType) {
  case NMD_NONE: break;
  case NMD_NMD: {
    Essex::CompositeNode *stopNode=
      new Essex::CompositeNode("premature-stop");
    stopNode->append("NMD");
    stopNode->append("EJC-distance",ejcDistance);
    status->append(stopNode);
  } break;
  //case NMD_TRUNCATION:  { // ### this is handled below
  case NMD_NO_STOP: 
    if(refTrans->hasUTR3()) status->append("nonstop-decay");
    break;
  case NMD_NO_START: status->append("no-start-codon"); break;
  }

  // Check for start codon
  if(startCodonMsg) status->append(startCodonMsg);

  // Check for frameshifts and amino acid differences
  const bool proteinsDiffer=refProtein!=altProtein;
  if(proteinsDiffer) {
    int matches;
    alignProteins(refProtein,altProtein,matches);
    Essex::CompositeNode *fate=new Essex::CompositeNode("protein-differs");
    percentMatch(matches,refProtein.length(),altProtein.length(),fate);
    status->append(fate);
    Labeling altLab(altSeqLen);
    computeLabeling(*altTrans,altLab);
    checker.checkFrameshifts(projectedLab,altLab,status);
  }

  // Check for premature stop codon
  if(nmdType==NMD_NONE && proteinsDiffer && refProtein.findFirst('*')>0 &&
     altProtein.findFirst('*')>0) {
    int refStop=refTrans->stopCodonGlobalCoord();
    int altStop=altTrans->stopCodonGlobalCoord();
    int altStopOnRef=(*revAlignment).mapApproximate(altStop,DIR_LEFT);
    if(altStopOnRef>=0 && refStop>=0 && altStopOnRef<refStop) {
      Essex::CompositeNode *fate=new Essex::CompositeNode("premature-stop");
      Essex::CompositeNode *truncation=
	new Essex::CompositeNode("protein-truncation");
      const int diff=getTruncationLength(*refTrans,altStopOnRef,refStop);
      truncation->append(String("")+diff+"aa");
      fate->append(truncation);
      status->append(fate);
    }
  }
  Essex::CompositeNode *uORFnode=
    orfAnalyzer->lostUORFs(*refTrans,refSeqStr,refSeq,*altTrans,altSeqStr,
			   altSeq,*alignment,reverseStrand);
  if(uORFnode) status->append(uORFnode);
}



/****************************************************************
 ACE::getTruncationLength()
 ****************************************************************/
int ACE::getTruncationLength(const GffTranscript &transcript,int PTC,int stop)
{
  // This works for the forward strand only!

  int nuc=0; bool accumulating=false;
  for(Vector<GffExon*>::const_iterator cur=transcript.getExons(), end=
	transcript.getExonsEnd() ; cur!=end ; ++cur) {
    GffExon &exon=**cur;
    if(exon.contains(PTC)) {
      accumulating=true;
      if(exon.contains(stop)) { nuc=stop-PTC; break; }
      else nuc=exon.getEnd()-PTC;
    }
    else if(accumulating) {
      if(exon.contains(stop)) { nuc+=stop-exon.getBegin(); break; }
      else nuc+=exon.getEnd()-exon.getBegin();
    }
  }
  return (nuc+2)/3;
}



/****************************************************************
 ACE::handleNoncoding()
 ****************************************************************/
void ACE::handleNoncoding(const GffTranscript *altTrans)
{
  if(!quiet) status->append("noncoding");
  int refOrfLen, altOrfLen; double refStartScore, altStartScore;
  String reason;
  Essex::CompositeNode *codingTranscript=
    orfAnalyzer->noncodingToCoding(*refTrans,refSeqStr,refSeq,*altTrans,
				   altSeqStr,altSeq,refOrfLen,altOrfLen,
				   reverseStrand,altSeqLen,refStartScore,
				   altStartScore,reason);
  if(codingTranscript) {
    Essex::CompositeNode *changeNode=
      new Essex::CompositeNode("noncoding-to-coding");
    changeNode->append("reason",reason);
    Essex::CompositeNode *lengthNode=
      new Essex::CompositeNode("ORF-length");
    lengthNode->append(refOrfLen);
    lengthNode->append("=>");
    lengthNode->append(altOrfLen);
    changeNode->append(lengthNode);
    changeNode->append("ref-start-score",float(refStartScore));
    changeNode->append("alt-start-score",float(altStartScore));
    const double cutoff=sensors.startCodonSensor->getCutoff();
    changeNode->append("start-score-cutoff",float(cutoff));
    changeNode->append(codingTranscript);
    status->append(changeNode);
  }
}



/****************************************************************
 ACE::enumerateAlts()
 ****************************************************************/
void ACE::enumerateAlts(Essex::CompositeNode *altTransEssex,
			TranscriptSignals *signals,
			GffTranscript *altTrans,
			ostream &osACE,
			const Labeling &projectedLab)
{
  altTransEssex->deleteChild("translation");
  appendBrokenSignals(signals);
  EnumerateAltStructures enumerator(*signals,altSeqStr,MAX_SPLICE_SHIFT,
				    MIN_EXON_LEN,MIN_INTRON_LEN,
				    NMD_DISTANCE_PARM,sensors,
				    allowExonSkipping,allowIntronRetention,
				    allowCrypticSites,reverseStrand);
  Vector<AlternativeStructure*> &altStructures=
    enumerator.getAltStructures();
  const int numStruct=altStructures.size();
  if(numStruct>0) {
    status->prepend("splicing-changes");
    Essex::CompositeNode *altStructNode=
      new Essex::CompositeNode("alternate-structures");
    status->append(altStructNode);
    for(Vector<AlternativeStructure*>::iterator cur=
	  altStructures.begin(), end=altStructures.end() ; cur!=end ; 
	++cur) {
      AlternativeStructure &s=**cur;
      processAltStructure(s,altStructNode,projectedLab);
    }
  }
  else status->prepend("no-transcript");
  delete altTrans;
}



/****************************************************************
 ACE::processAltStructure()
 ****************************************************************/
Essex::CompositeNode *ACE::processAltStructure(AlternativeStructure &s,
			      Essex::CompositeNode *altStructNode,
			      const Labeling &projectedLab)
{
  Essex::CompositeNode *msg2=NULL, *changeToCoding=NULL, *uORFnode=NULL;
  if(refTrans->isCoding()) { // see if a new start codon extends the ORF
    int oldOrfLen, newOrfLen; double oldStartScore, newStartScore;
    String oldStartStr, newStartStr; Essex::CompositeNode *reason;
    GffTranscript *newTranscript=
      orfAnalyzer->earlierStartCodon(*refTrans,refSeqStr,refSeq,
				     *s.transcript,altSeqStr,altSeq,
				     *revAlignment,oldOrfLen,newOrfLen,
				     oldStartScore,newStartScore,
				     oldStartStr, newStartStr,
				     reverseStrand,altSeqLen,reason);
    if(newTranscript) {
      if(reverseStrand) { // added 7/28/2016
	GffTranscript temp(*newTranscript);
	temp.reverseComplement(altSeqLen);
	changeToCoding=temp.toEssex();
      } else changeToCoding=newTranscript->toEssex();
      msg2=
	new Essex::CompositeNode("new-upstream-start-codon");
      const double cutoff=sensors.startCodonSensor->getCutoff();
      msg2->append(signalNode("new-start",newStartStr,newStartScore,cutoff));
      msg2->append(signalNode("old-start",oldStartStr,oldStartScore,cutoff));
      msg2->append(reason);
      Essex::CompositeNode *lengthNode=
	new Essex::CompositeNode("ORF-length");
      lengthNode->append(oldOrfLen);
      lengthNode->append("=>");
      lengthNode->append(newOrfLen);
      msg2->append(lengthNode);
      Essex::CompositeNode *fate=new Essex::CompositeNode("fate");
      ProjectionChecker checker(*refTrans,*newTranscript,refSeqStr,refSeq,
				altSeqStr,altSeq,projectedLab,sensors);
      analyzeEarlierStart(newTranscript,checker,projectedLab,fate);
      changeToCoding->append(fate);
      delete s.transcript;
      s.transcript=newTranscript;
    }
    uORFnode=orfAnalyzer->lostUORFs(*refTrans,refSeqStr,refSeq,*s.transcript,
				    altSeqStr,altSeq,*alignment,reverseStrand);
  }
  else { // noncoding: check whether it changes to coding
    int refOrfLen, altOrfLen;
    double refStartScore, altStartScore; String reason;
    changeToCoding=
      orfAnalyzer->noncodingToCoding(*refTrans,refSeqStr,refSeq,*s.transcript,
				     altSeqStr,altSeq,refOrfLen,altOrfLen,
				     reverseStrand,altSeqLen,refStartScore,
				     altStartScore,reason);
    if(changeToCoding) {
      msg2=new Essex::CompositeNode("noncoding-to-coding");
      msg2->append("reason",reason);
      Essex::CompositeNode *lengthNode=
	new Essex::CompositeNode("ORF-length");
      lengthNode->append(refOrfLen);
      lengthNode->append("=>");
      lengthNode->append(altOrfLen);
      msg2->append(lengthNode);
      msg2->append("ref-start-score",float(refStartScore));
      msg2->append("alt-start-score",float(altStartScore));
      const double cutoff=sensors.startCodonSensor->getCutoff();
      msg2->append("start-score-cutoff",float(cutoff));
    }
  }
  Essex::CompositeNode *node=NULL;
  if(changeToCoding) node=changeToCoding;
  else {
    s.transcript->loadSequence(altSeqStr);
    s.transcript->computePhases();
    GffTranscript temp(*s.transcript);
    if(reverseStrand) temp.reverseComplement(altSeqLen);
    node=temp.toEssex();
  }
  if(uORFnode) node->append(uORFnode);
  VariantClassifier classifier(variants,VariantClassifier::ALT,
			       *s.transcript,altSeqLen,reverseStrand);
  node->append(classifier.makeVariantsNode());
  s.reportCrypticSites(node,reverseStrand,altSeqLen);
  listStructureChanges(s,node,s.msg,msg2);
  if(s.msg) { status->append(s.msg); s.msg=NULL; }
  altStructNode->append(node);
  if(!changeToCoding) handleProteinFate(s,node);
  return node;
}



/****************************************************************
 ACE::listStructureChanges()
 ****************************************************************/
void ACE::listStructureChanges(const AlternativeStructure &s,
			       Essex::CompositeNode *node,
			       Essex::CompositeNode *&msg,
			       Essex::CompositeNode *msg2)
{
  if(s.structureChange.anyChange()) {
    Essex::CompositeNode *changeNode=
      new Essex::CompositeNode("structure-change");
    node->prepend(changeNode);
    if(s.structureChange.crypticSite) 
      changeNode->append("cryptic-site");
    if(s.structureChange.exonSkipping) 
      changeNode->append("exon-skipping");
    if(s.structureChange.intronRetention) 
      changeNode->append("intron-retention");
    if(s.structureChange.deNovoSite) 
      changeNode->append("de-novo-site");
    if(s.structureChange.crypticExon) 
      changeNode->append("cryptic-exon");
    if(s.structureChange.regulatoryChange) 
      changeNode->append("splicing-regulatory-change");
    if(msg) { changeNode->append(msg); msg=NULL; }
    if(msg2) changeNode->append(msg2);
  }
  //else throw "no changes";
  else {
    Essex::CompositeNode *changeNode=
      new Essex::CompositeNode("structure-change");
    changeNode->append("mapped-transcript");
    node->prepend(changeNode);
  }
}



/****************************************************************
 ACE::handleProteinFate()
 ****************************************************************/
void ACE::handleProteinFate(const AlternativeStructure &s,
			    Essex::CompositeNode *node)
{
  switch(s.proteinFate) {
  case NMD_NONE:{      // nothing wrong
    s.transcript->loadSequence(altSeqStr);
    bool identical=s.transcript->getProtein()==refProtein;
    if(identical)
      node->append("fate","identical-protein");
    else {
      int matches, len;
      String altProtein=s.transcript->getProtein();
      alignProteins(refProtein,altProtein,matches);
      Essex::CompositeNode *fate=new Essex::CompositeNode("fate");
      fate->append("protein-differs");
      percentMatch(matches,refProtein.length(),altProtein.length(),
		   fate);
      node->append(fate);
    }
  }break;
  case NMD_NMD: {       // premature stop codon & NMD
    Essex::CompositeNode *fate=new Essex::CompositeNode("fate");
    fate->append("NMD");
    fate->append("EJC-distance",s.ejcDistance); // ### always zero -- why?
    node->append(fate);
  }break;
  case NMD_TRUNCATION: {// premature stop codon, truncated protein
    int matches, len;
    String altProtein=s.transcript->getProtein();
    alignProteins(refProtein,altProtein,matches);
    Essex::CompositeNode *fate=new Essex::CompositeNode("fate");
    fate->append("protein-truncation");
    percentMatch(matches,refProtein.length(),altProtein.length(),fate);
    node->append(fate);
  }
    break;
  case NMD_NO_STOP:    // no stop codon
    if(refTrans->hasUTR3()) // can't predict if no annotated UTR
      node->append("fate","nonstop-decay");
    break;
  case NMD_NO_START:   // no start codon
    node->append("fate","noncoding");
    break;
  }
}



/****************************************************************
 ACE::appendBrokenSignals()
 ****************************************************************/
void ACE::appendBrokenSignals(const TranscriptSignals *signals)
{
  int numSignals=signals->numSignals();
  for(int i=0 ; i<numSignals ; ++i) {
    TranscriptSignal &signal=(*signals)[i];
    if(!signal.isBroken() && !signal.isWeakened()) continue;
    String tag;
    SignalType type=signal.getType();
    if(type==GT) tag=signal.weakened ? "weakened-donor" : "broken-donor";
    else if(type==AG) 
      tag=signal.weakened ? "weakened-acceptor" : "broken-acceptor";
    else INTERNAL_ERROR;
    Essex::CompositeNode *node=new Essex::CompositeNode(tag);
    int pos=signal.getPos();
    if(reverseStrand) pos=altSeqLen-pos-2; // -2 is for the signal length!
    node->append(pos);
    Vector<String> fields; signal.seq.getFields(fields);
    for(Vector<String>::iterator cur=fields.begin(), end=fields.end() ; 
	cur!=end ; ++cur) node->append(*cur);
    SignalSensor *sensor=sensors.findSensor(type);
    node->append("threshold:"); node->append(float(sensor->getCutoff()));
    status->append(node);
  }
}



/****************************************************************
 ACE::writeXML()
 ****************************************************************/
void ACE::writeXML()
{
  ofstream os(xmlFilename.c_str());
  root->printXML(os);
  os<<endl;
}



/****************************************************************
 ACE::append()
 ****************************************************************/
void ACE::append(Essex::CompositeNode *root,const String &tag,
			 const String &message)
{
  Essex::CompositeNode *node=new Essex::CompositeNode(tag);
  node->append(message);
  root->appendChild(node);
}



/****************************************************************
 ACE::append()
 ****************************************************************/
void ACE::append(Essex::CompositeNode *root,const char *tag,
			 const char *message)
{
  append(root,String(tag),String(message));
}



/****************************************************************
 ACE::append()
 ****************************************************************/
void ACE::append(Essex::CompositeNode *root,const char *tag,int x)
{
  Essex::CompositeNode *node=new Essex::CompositeNode(tag);
  node->append(x);
  root->appendChild(node);
}



/****************************************************************
 ACE::append()
 ****************************************************************/
void ACE::append(Essex::CompositeNode *root,const char *tag,
			const String &message)
{
  append(root,tag,message.c_str());
}



/****************************************************************
 ACE::computeLabeling()
 ****************************************************************/
void ACE::computeLabeling(GffTranscript &transcript,
				  Labeling &refLab)
{
  const int begin=transcript.getBegin(), end=transcript.getEnd();
  char strand=transcript.getStrand();
  if(strand!='+') throw "only forward-strand features are currently supported";
  int numExons=transcript.getNumExons();
  refLab.asArray().setAllTo(LABEL_INTERGENIC);
  for(int i=begin ; i<end ; ++i) refLab[i]=LABEL_INTRON;
  int phase=0;
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript.getIthExon(i);
    const int begin=exon.getBegin(), end=exon.getEnd();
    for(int j=begin ; j<end ; ++j) {
      refLab[j]=getExonLabel(phase);
      phase=(phase+1)%3;
    }
  }
  for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	transcript.getUTRend() ; cur!=end ; ++cur) {
    GffExon *UTR=*cur;
    refLab.setIntervalTo(Interval(UTR->getBegin(),UTR->getEnd()),LABEL_UTR);
  }
}



/****************************************************************
 ACE::mapLabeling()
 ****************************************************************/
void ACE::mapLabeling(Labeling &from,Labeling &to,const String &filename)
{
  CigarAlignment &align=*alignment;
  to.asArray().setAllTo(LABEL_NONE);
  int L=align.length();
  for(int i=0 ; i<L ; ++i) {
    int j=align[i];
    if(j!=CIGAR_UNDEFINED) to[j]=from[i];
  }
  if(!filename.empty()) {
    ofstream os(filename.c_str());
    os<<to;
    os.close();
  }
}



/****************************************************************
 ACE::mapExon()
 ****************************************************************/
bool ACE::mapExon(GffExon &exon,CigarAlignment &align)
{
  int begin=exon.getBegin(), end=exon.getEnd();

  // These two lines map the splice sites across the
  // alignment, then use that to set exon boundaries:
  begin=align.mapApproximate(begin-2,DIR_NONE)+2;
  end=align.mapApproximate(end,DIR_NONE);
  if(begin<0 || end<0) return false;
  exon.setBegin(begin); exon.setEnd(end);
  return true;
}



/****************************************************************
 ACE::mapTranscript()
 ****************************************************************/
void ACE::mapTranscript(const String &outfile)
{
  CigarAlignment &align=*alignment;
  Vector<GffExon*> rawExons;
  refTrans->getRawExons(rawExons);
  GffTranscript transcript(refTrans->getTranscriptId(),
			   refTrans->getSubstrate(),
			   refTrans->getStrand(),"ACE");
  transcript.setGeneId(refTrans->getGeneId());
  transcript.setSubstrate(substrate);
  transcript.getSource()="ACE";
  for(Vector<GffExon*>::iterator cur=rawExons.begin(), end=rawExons.end() ;
	cur!=end ; ++cur) {
    GffExon *exon=new GffExon(**cur,transcript);
    if(!mapExon(*exon,align)) INTERNAL_ERROR;
    transcript.addUTR(exon);
  }
  GffTranscript::deleteExons(rawExons);
  if(refTrans->isCoding()) {
    int mappedStartCodon=
      align.mapApproximate(refTrans->getIthExon(0).getBegin(),DIR_LEFT);
    int startCodon=StartCodonFinder::findStartCodon(transcript,
						    transcript.peekUTR(),
						    altSeqStr,
						    mappedStartCodon,
						    sensors);
    if(startCodon!=mappedStartCodon)
      if(startCodon>0) {
	Essex::CompositeNode *newNode
	  =new Essex::CompositeNode("start-codon-change");
	if(reverseStrand) {
	  newNode->append("from",altSeqLen-mappedStartCodon);
	  newNode->append("to",altSeqLen-startCodon); }
	else {
	  newNode->append("from",mappedStartCodon);
	  newNode->append("to",startCodon); }
	startCodonMsg=newNode; }
      else startCodonMsg=new Essex::StringNode("no-start-codon");
    if(startCodon>=0)
      transcript.splitUTRandCDS(altSeqStr,startCodon,sensors.stopCodons);
  }
  transcript.setExonTypes(); transcript.setUTRtypes();
  ofstream os(outfile.c_str());
  transcript.toGff(os);
}



/****************************************************************
 ACE::loadSeq()
 ****************************************************************/
String ACE::loadSeq(const String &filename)
{
  FastaReader reader(filename);
  String def, seq;
  if(!reader.nextSequence(def,seq)) throw filename+" : cannot read file";
  return seq;
}



/****************************************************************
 ACE::loadSeq()
 ****************************************************************/
String ACE::loadSeq(const String &filename,String &CIGAR)
{
  FastaReader reader(filename);
  String seq, remainder;
  if(!reader.nextSequence(altDefline,seq)) 
    throw filename+" : cannot read file";
  const int L=seq.length();
  FastaReader::parseDefline(altDefline,substrate,remainder);
  if(warningsRegex.search(remainder)) VCFwarnings=warningsRegex[1];
  if(errorsRegex.search(remainder)) VCFerrors=errorsRegex[1];
  if(coordRegex.search(remainder)) globalCoord=coordRegex[1];
  Map<String,String> attr;
  FastaReader::parseAttributes(remainder,attr);
  if(!attr.isDefined("cigar")) 
    throw String("No CIGAR string found on defline: ")+altDefline;
  CIGAR=attr["cigar"];
  parseVariants(attr["variants"],variants,L);
  return seq;
}



/****************************************************************
 ACE::loadGff()
 ****************************************************************/
GffTranscript *ACE::loadGff(const String &filename)
{
  GffReader reader(filename);
  Vector<GffTranscript*> *transcripts=reader.loadTranscripts();
  const int n=transcripts->size();
  if(n<1) throw filename+" contains no transcripts";
  GffTranscript *transcript=(*transcripts)[0];
  for(int i=1 ; i<n ; ++i) delete (*transcripts)[i];
  delete transcripts;
  transcript->setExonTypes();
  transcript->setUTRtypes();
  return transcript;
}



/****************************************************************
 ACE::writeProtein()
 ****************************************************************/
void ACE::writeProtein(const String &def,const String &protein,
			       const String &filename)
{
  if(filename.empty()) return;
  fastaWriter.writeFasta(def,protein,filename);
}



/****************************************************************
 ACE::processConfig()
 ****************************************************************/
void ACE::processConfig(const String &filename)
{
  const String path=File::getPath(filename);
  char *oldPath=new char[PATH_MAX];
  getcwd(oldPath,PATH_MAX);
  chdir(path.c_str());

  ConfigFile config(filename);
  MAX_SPLICE_SHIFT=config.getIntOrDie("max-splice-shift");
  MIN_EXON_LEN=config.getIntOrDie("min-exon-length");
  MIN_INTRON_LEN=config.getIntOrDie("min-intron-length");
  NMD_DISTANCE_PARM=config.isDefined("NMD-distance") ?
    config.getIntOrDie("NMD-distance") : 50;
  nmd.setDistParm(NMD_DISTANCE_PARM);

  openPenalty=-config.getFloatOrDie("gap-open-penalty");
  extendPenalty=-config.getFloatOrDie("gap-extend-penalty");
  bandwidth=config.getIntOrDie("bandwidth");
  String matrixFile=config.lookupOrDie("subst-matrix");
  substMatrix=
    new SubstitutionMatrix<float>(matrixFile,AminoAlphabet::global());

  allowExonSkipping=config.getBoolOrDie("allow-exon-skipping");
  allowIntronRetention=config.getBoolOrDie("allow-intron-retention");
  allowCrypticSites=config.getBoolOrDie("allow-cryptic-sites");

  parseConsensusList("donor-consensus",config,sensors.donorConsensuses);
  parseConsensusList("acceptor-consensus",config,sensors.acceptorConsensuses);
  parseConsensusList("stop-codons",config,sensors.stopCodons);
  parseConsensusList("start-codons",config,sensors.startCodons);

  sensors.startCodonSensor=loadModel("start-codon-model",config);
  sensors.shortStartSensor=loadModel("short-start-codon-model",config);
  sensors.stopCodonSensor=loadModel("stop-codon-model",config);
  sensors.donorSensor=loadModel("donor-model",config);
  sensors.acceptorSensor=loadModel("acceptor-model",config);
  sensors.setConsensuses();

  const int MIN_ORF_LEN=config.getIntOrDie("min-orf-length");
  orfAnalyzer=new OrfAnalyzer(sensors,MIN_ORF_LEN);

  chdir(oldPath);
  delete [] oldPath;
}



/****************************************************************
 ACE::parseConsensusList()
 ****************************************************************/
void ACE::parseConsensusList(const String &tag,ConfigFile &config,
				     Set<String> &into)
{
  String consensusString=config.lookupOrDie(tag);
  Vector<String> fields;
  consensusString.getFields(fields,",");
  for(Vector<String>::const_iterator cur=fields.begin(), end=fields.end() ;
      cur!=end ; ++cur)
    into.insert(*cur);
}



/****************************************************************
 ACE::loadModel()
 ****************************************************************/
SignalSensor *ACE::loadModel(const String &label,ConfigFile &config)
{
  String filename=config.lookupOrDie(label);
  return SignalSensor::load(filename,garbageCollector);
}



/****************************************************************
 ACE::alignProteins()
 ****************************************************************/
double ACE::alignProteins(const String &refStr,const String &altStr,
			 int &matches)
{
  if(refStr.length()==0 || altStr.length()==0) return 0.0;
  const AminoAlphabet &alphabet=AminoAlphabet::global();
  Sequence refSeq(refStr,alphabet), altSeq(altStr,alphabet);
  BandedSmithWaterman<float> aligner(alphabet,refSeq,altSeq,*substMatrix,
				     openPenalty,extendPenalty,bandwidth);
  Alignment *alignment=aligner.fullAlignment();
  matches=alignment->countMatches();
  double score=double(matches)/refStr.length();
  delete alignment;
  return score;
}



/****************************************************************
 ACE::percentMatch()
 ****************************************************************/
void ACE::percentMatch(int matches,int refLen,int altLen,
		       Essex::CompositeNode *parent)
{
  double percent=double(matches)/max(refLen,altLen);
  Essex::CompositeNode *node=new Essex::CompositeNode("percent-match");
  node->append(String(int(10000*percent+5.0/9)/100.0));
  node->append(String(matches)+"/"+max(refLen,altLen));
  node->append(String("ref-length=")+refLen);
  node->append(String("alt-length=")+altLen);
  parent->append(node);
}



/****************************************************************
 ACE::parseVariants()
 ****************************************************************/
void ACE::parseVariants(const String &s,Vector<Variant> &variants,int L)
{
  Vector<String> fields;
  s.getFields(fields,",");
  const int numFields=fields.size();
  for(int i=0 ; i<numFields ; ++i) {
    if(!variantRegex.match(fields[i])) throw "Can't parse variant "+fields[i];
    String id=variantRegex[1], chr=variantRegex[2];
    int refPos=variantRegex[3].asInt(), altPos=variantRegex[4].asInt();
    String ref=variantRegex[5], alt=variantRegex[6];
    Variant v(id,chr,refPos,altPos,i);
    v.addAllele(ref); v.addAllele(alt);
    variants.push_back(v);
  }
}



/****************************************************************
 ACE::makeEssexVariants()
 ****************************************************************/
Essex::CompositeNode *ACE::makeEssexVariants()
{
  Essex::CompositeNode *parent=new Essex::CompositeNode("variants");
  const int numVariants=variants.size();
  for(int i=0 ; i<numVariants ; ++i) {
    const Variant &v=variants[i];
    String s=v.id+":"+v.chr+":"+v.refPos+":"+v.altPos+":"+v.alleles[0]
      +":"+v.alleles[1];
    parent->append(s);
  }
  return parent;
}



