/****************************************************************
 ACE.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
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
#include "Labeling.H"
#include "ProjectionChecker.H"
#include "EnumerateAltStructures.H"
#include "SignalSensors.H"
#include "GarbageCollector.H"
#include "StartCodonFinder.H"
using namespace std;
using namespace BOOM;





ACE::ACE()
  : warningsRegex("/warnings=(\\d+)"), errorsRegex("/errors=(\\d+)"), 
    VCFwarnings(0), VCFerrors(0), startCodonMsg(NULL)
  {
    // ctor
  }



ACE::~ACE()
{
  cout<<"ACE terminated successfully"<<endl;
}



int ACE::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"a:cd:e:l:x:");
  if(cmd.numArgs()!=6)
    throw String("\n\
ace <ace.config> <ref.gff> <ref.fasta> <alt.fasta> <out.gff> <out.essex>\n\
     -a acceptors = also allow these acceptors; comma-separated (capital)\n\
     -c = sequence has been reversed, but cigar string has not\n\
     -d donors = also allow these donors; comma-separated (capital letters)\n\
     -e N = abort if vcf errors >N\n\
     -l <file> = emit a per-nucleotide labeling for the alt sequence\n\
     -x <file> = also emit xml\n\
  alt.fasta must have a cigar string: >ID ... /cigar=1045M3I10M7D4023M ...\n\
\n");
  const String configFile=cmd.arg(0);
  const String refGffFile=cmd.arg(1);
  const String refFasta=cmd.arg(2);
  const String altFasta=cmd.arg(3);
  const String outGff=cmd.arg(4);
  const String outACE=cmd.arg(5);
  alphabet=DnaAlphabet::global();
  if(cmd.option('l')) labelingFile=cmd.optParm('l');
  if(cmd.option('x')) xmlFilename=cmd.optParm('x');
  reverseCigar=cmd.option('c');

  // Read some data from files
  processConfig(configFile);
  String refSeqStr=loadSeq(refFasta), altSeqStr=loadSeq(altFasta,CIGAR);
  const Sequence refSeq(refSeqStr,alphabet), altSeq(altSeqStr,alphabet);
  int refSeqLen=refSeqStr.length(), altSeqLen=altSeqStr.length();
  GffTranscript *refTrans=loadGff(refGffFile);
  refTrans->loadSequence(refSeqStr);
  String refProtein=refTrans->getProtein();

  // Set up to generate structured output in Essex/XML
  String transcriptID=refTrans->getTranscriptId();
  String geneID=refTrans->getGeneId();
  root=new Essex::CompositeNode("report");
  append(root,"substrate",substrate);
  append(root,"transcript-ID",transcriptID);
  append(root,"gene-ID",geneID);
  append(root,"vcf-warnings",VCFwarnings);
  append(root,"vcf-errors",VCFerrors);
  append(root,"alignment",CIGAR);
  append(root,"defline",altDefline);
  Essex::CompositeNode *refTransEssex=refTrans->toEssex();
  refTransEssex->getTag()="reference-transcript";
  root->append(refTransEssex);
  Essex::CompositeNode *status=new Essex::CompositeNode("status");
  root->appendChild(status);
  ofstream osACE(outACE.c_str());
  if(cmd.option('e') && VCFerrors>cmd.optParm('e').asInt()) {
    status->append("too-many-vcf-errors");
    if(!xmlFilename.empty()) writeXML();
    osACE<<*root<<endl;
    osACE<<"#===========================================================\n";
    return -1;
  }
  
  // Check that the reference gene is well-formed
  bool noStart, noStop, PTC, badSpliceSite, referenceIsOK=true;
  String msg;
  if(!ProjectionChecker::geneIsWellFormed(*refTrans,refSeqStr,
					  noStart,noStop,PTC,badSpliceSite,
					  status,sensors)) {
    status->prepend("bad-annotation");
    referenceIsOK=false;
  }

  // Compute the reference labeling
  Labeling refLab(refSeqLen);
  computeLabeling(*refTrans,refLab);

  // Project the reference labeling over to the alternate sequence
  Labeling altLab(altSeqLen);
  CigarString cigar(CIGAR);
  if(reverseCigar) cigar.reverse();
  mapLabeling(refLab,altLab,cigar);

  // Project the reference GFF over to an alternate GFF
  mapTranscript(*refTrans,cigar,outGff,altSeqStr,altSeq);

  // Generate labeling file
  if(!labelingFile.empty()) {
    ofstream os(labelingFile.c_str());
    os<<altLab;
    os.close(); }

  // Check the projection to see if the gene might be broken
  if(referenceIsOK) {
    GffTranscript *altTrans=loadGff(outGff);
    altTrans->loadSequence(altSeqStr);
    Essex::CompositeNode *altTransEssex=altTrans->toEssex();
    altTransEssex->getTag()="mapped-transcript";
    root->append(altTransEssex);
    ProjectionChecker checker(*refTrans,*altTrans,refSeqStr,refSeq,
			      altSeqStr,altSeq,altLab,sensors);
    TranscriptSignals *signals=checker.findBrokenSpliceSites();
    if(!signals) {
      status->prepend("unequal-numbers-of-exons");
      if(!xmlFilename.empty()) writeXML();
      osACE<<*root<<endl;
      osACE<<"#===========================================================\n";
      return -1;
    }

    if(signals->anyBroken()) {
      appendBrokenSignals(signals,status);
      EnumerateAltStructures enumerator(*signals,altSeqStr,MAX_SPLICE_SHIFT,
					MIN_EXON_LEN,MIN_INTRON_LEN,sensors,
					allowExonSkipping,allowIntronRetention,
					allowCrypticSites);
      const Vector<AlternativeStructure*> &altStructures=
	enumerator.getAltStructures();
      const int numStruct=altStructures.size();
      if(numStruct>0) {
	status->prepend("splicing-changes");
	Essex::CompositeNode *altStructNode=
	  new Essex::CompositeNode("alternate-structures");
	status->append(altStructNode);
	for(Vector<AlternativeStructure*>::const_iterator cur=
	      altStructures.begin(), end=altStructures.end() ; cur!=end ; 
	    ++cur) {
	  const AlternativeStructure &s=**cur;
	  Essex::CompositeNode *msg=s.msg;
	  s.transcript->loadSequence(altSeqStr);
	  Essex::CompositeNode *node=s.transcript->toEssex();
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
	    if(msg) { changeNode->append(msg); msg=NULL; }
	  }
	  if(msg) status->append(msg);
	  altStructNode->append(node);
	  switch(s.proteinFate) {
	  case NMD_NONE:{      // nothing wrong
	    s.transcript->loadSequence(altSeqStr);
	    bool identical=s.transcript->getProtein()==refProtein;
	    node->append("fate",
			 identical ? "identical-protein" : "protein-differs");
	  }break;
	  case NMD_NMD:        // premature stop codon & NMD
	    node->append("fate","NMD");
	    break;
	  case NMD_TRUNCATION: // premature stop codon, truncated protein
	    node->append("fate","protein-truncation");
	    break;
	  case NMD_NO_STOP:    // no stop codon
	    node->append("fate","nonstop-decay");
	    break;
	  case NMD_NO_START:   // no start codon
	    node->append("fate","noncoding");
	    break;
	  }
	  
	}
      }
      else status->prepend("no-transcript");
      if(!xmlFilename.empty()) writeXML();
      osACE<<*root<<endl;
      osACE<<"#===========================================================\n";
      return 0;
    }

    // Otherwise, projection was successful
    status->prepend("mapped");

    // Translate to proteins
    if(refTrans->isCoding()) {
      String refProtein, altProtein;
      checker.translate(*refTrans,*altTrans,refProtein,altProtein);
      if(!refProteinFile.empty())
	writeProtein(">ref",refProtein,refProteinFile);
      if(!altProteinFile.empty())
	writeProtein(">alt",altProtein,altProteinFile);
      
      // Check for start codon
      if(startCodonMsg) status->append(startCodonMsg);
      else if(!checker.hasStartCodon(altProtein)) 
	status->append("no-start-codon");
      
      // Check for frameshifts
      if(refProtein==altProtein) status->append("identical-protein");
      else {
	status->append("protein-differs");
	checker.checkFrameshifts(altLab,*altTrans,altSeqStr,status); }
      
      // Check for stop codons
      int PTCpos;
      if(checker.hasPTC(altProtein,PTCpos)) {
	Essex::CompositeNode *node=new Essex::CompositeNode("premature-stop");
	int EJCdistance;
	bool nmd=checker.detectNMD(*altTrans,altSeqStr,false,EJCdistance);
	node->append(nmd ? "NMD" : "protein-truncation");
	node->append("EJC-distance",EJCdistance);
	node->append("AA-pos",PTCpos);
	status->append(node);
      }
      else if(!checker.hasStopCodon(altProtein)) {
	refTrans->extendFinalExonBy3(); altTrans->extendFinalExonBy3();
	checker.translate(*refTrans,*altTrans,refProtein,altProtein);
	if(!checker.hasStopCodon(altProtein)) status->append("nonstop-decay");
      }
    }
    else status->append("noncoding");
  }

  // Flush output
  if(!xmlFilename.empty()) writeXML();
  osACE<<*root<<endl;
  osACE<<"#===========================================================\n";
  return 0;
}



void ACE::appendBrokenSignals(const TranscriptSignals *signals,
				      Essex::CompositeNode *status)
{
  int numSignals=signals->numSignals();
  for(int i=0 ; i<numSignals ; ++i) {
    TranscriptSignal &signal=(*signals)[i];
    if(!signal.isBroken()) continue;
    String tag;
    SignalType type=signal.getType();
    if(type==GT) tag=signal.weakened ? "weakened-donor" : "broken-donor";
    else if(type==AG) 
      tag=signal.weakened ? "weakened-acceptor" : "broken-acceptor";
    else INTERNAL_ERROR;
    Essex::CompositeNode *node=new Essex::CompositeNode(tag);
    node->append(signal.getPos());
    Vector<String> fields; signal.seq.getFields(fields);
    for(Vector<String>::iterator cur=fields.begin(), end=fields.end() ; 
	cur!=end ; ++cur) node->append(*cur);
    status->append(node);
  }
}



void ACE::writeXML()
{
  ofstream os(xmlFilename.c_str());
  root->printXML(os);
  os<<endl;
}



void ACE::append(Essex::CompositeNode *root,const String &tag,
			 const String &message)
{
  Essex::CompositeNode *node=new Essex::CompositeNode(tag);
  node->append(message);
  root->appendChild(node);
}



void ACE::append(Essex::CompositeNode *root,const char *tag,
			 const char *message)
{
  append(root,String(tag),String(message));
}



void ACE::append(Essex::CompositeNode *root,const char *tag,int x)
{
  Essex::CompositeNode *node=new Essex::CompositeNode(tag);
  node->append(x);
  root->appendChild(node);
}



void ACE::append(Essex::CompositeNode *root,const char *tag,
			const String &message)
{
  append(root,tag,message.c_str());
}



void ACE::computeLabeling(GffTranscript &transcript,
				  Labeling &refLab)
{
  //int begin, end;
  //transcript.getCDSbeginEnd(begin,end);
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



void ACE::mapLabeling(Labeling &from,Labeling &to,
			      const CigarString &cigar)
{
  CigarAlignment &align=*cigar.getAlignment();
  to.asArray().setAllTo(LABEL_NONE);
  int L=align.length();
  for(int i=0 ; i<L ; ++i) {
    int j=align[i];
    if(j!=CIGAR_UNDEFINED) to[j]=from[i];
  }
  delete &align;
}



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



/*
void ACE::mapTranscript(GffTranscript &refTrans,
				const CigarString &cigar,
				const String &outfile,
				const String &altSeqStr,
				const Sequence &altSeq)
{
  GffTranscript transcript=refTrans;
  transcript.setSubstrate(substrate);
  transcript.getSource()="ACE";
  for(Vector<GffExon*>::iterator cur=transcript.getExons(), end=
	transcript.getExonsEnd() ; cur!=end ; ++cur)
    mapExon(**cur,align);
  for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	transcript.getUTRend() ; cur!=end ; ++cur) {
    mapExon(**cur,align);
  }
  delete &align;
  ofstream os(outfile.c_str());
  transcript.toGff(os);
}
 */



void ACE::mapTranscript(GffTranscript &refTrans,
				const CigarString &cigar,
				const String &outfile,
				const String &altSeqStr,
				const Sequence &altSeq)
{
  CigarAlignment &align=*cigar.getAlignment();
  Vector<GffExon*> rawExons;
  refTrans.getRawExons(rawExons);
  GffTranscript transcript(refTrans.getTranscriptId(),
			   refTrans.getSubstrate(),
			   refTrans.getStrand(),"ACE");
  transcript.setGeneId(refTrans.getGeneId());
  transcript.setSubstrate(substrate);
  transcript.getSource()="ACE";
  for(Vector<GffExon*>::iterator cur=rawExons.begin(), end=rawExons.end() ;
	cur!=end ; ++cur) {
    GffExon *exon=new GffExon(**cur,transcript);
    if(!mapExon(*exon,align)) INTERNAL_ERROR;
    transcript.addUTR(exon);
  }
  GffTranscript::deleteExons(rawExons);
  if(refTrans.isCoding()) {
    int mappedStartCodon=
      align.mapApproximate(refTrans.getIthExon(0).getBegin(),DIR_LEFT);
    int startCodon=StartCodonFinder::findStartCodon(transcript,
						    transcript.peekUTR(),
						    altSeqStr,
						    mappedStartCodon,
						    sensors);
    if(startCodon!=mappedStartCodon)
      if(startCodon>0) {
	startCodonMsg=new Essex::CompositeNode("start-codon-change");
	startCodonMsg->append("from",mappedStartCodon);
	startCodonMsg->append("to",startCodon); }
      else startCodonMsg=new Essex::CompositeNode("start-codon-lost");
    if(startCodon>=0)
      transcript.splitUTRandCDS(altSeqStr,startCodon,sensors.stopCodons);
  }
  transcript.setExonTypes(); transcript.setUTRtypes();
  delete &align;
  ofstream os(outfile.c_str());
  transcript.toGff(os);
}



String ACE::loadSeq(const String &filename)
{
  FastaReader reader(filename);
  String def, seq;
  if(!reader.nextSequence(def,seq)) throw filename+" : cannot read file";
  return seq;
}



String ACE::loadSeq(const String &filename,String &CIGAR)
{
  FastaReader reader(filename);
  String seq, remainder;
  if(!reader.nextSequence(altDefline,seq)) 
    throw filename+" : cannot read file";
  FastaReader::parseDefline(altDefline,substrate,remainder);
  if(warningsRegex.search(remainder)) VCFwarnings=warningsRegex[1];
  if(errorsRegex.search(remainder)) VCFerrors=errorsRegex[1];
  Map<String,String> attr;
  FastaReader::parseAttributes(remainder,attr);
  if(!attr.isDefined("cigar")) 
    throw String("No CIGAR string found on defline: ")+altDefline;
  CIGAR=attr["cigar"];
  return seq;
}



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



void ACE::writeProtein(const String &def,const String &protein,
			       const String &filename)
{
  if(filename.empty()) return;
  fastaWriter.writeFasta(def,protein,filename);
}



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

  chdir(oldPath);
  delete [] oldPath;
}



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



SignalSensor *ACE::loadModel(const String &label,ConfigFile &config)
{
  String filename=config.lookupOrDie(label);
  return SignalSensor::load(filename,garbageCollector);
}




