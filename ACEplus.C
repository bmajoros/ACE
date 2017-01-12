/****************************************************************
 ACEplus.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ACEplus.H"
#include "GraphBuilder.H"
#include "TranscriptPaths.H"
#include "GeometricDistribution.H"
#include "EmpiricalDistribution.H"
using namespace std;
using namespace BOOM;

/****************************************************************
 Globals
 ****************************************************************/
static bool VERBOSE=true;
static const char *PROGRAM_NAME="ACEplus";
static const char *VERSION="1.0";



/****************************************************************
 ACEplus::ACEplus()
 ****************************************************************/
ACEplus::ACEplus()
{
  // ctor
}



/****************************************************************
 ACEplus::main()
 ****************************************************************/
int ACEplus::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"ce:l:qx:v");
  parseCommandLine(cmd);

  // Read some data from files
  if(VERBOSE) cerr<<"loading inputs"<<endl;
  loadInputs(configFile,refGffFile,refFasta,altFasta);

  // Build prefix-sum arrays for fast scoring
  buildPSAs();

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
    checkProjection(outGff,mapped,refLab,projectedLab,osACE);

  // Flush output
  if(VERBOSE) cerr<<"cleaning up"<<endl;
  flushOutput(osACE,mapped);

  cout<<"ACE terminated successfully"<<endl;
  return 0;
}



/****************************************************************
 ACEplus::parseCommandLine()
 ****************************************************************/
void ACEplus::parseCommandLine(const CommandLine &cmd)
{
  if(cmd.option('v')) { cout<<"version "<<VERSION<<endl; exit(0); }
  if(cmd.numArgs()!=6)
    throw String("\n\
aceplus <ace.config> <ref.gff> <ref.fasta> <alt.fasta> <out.gff> <out.essex>\n\
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
 ACEplus::processConfig()
 ****************************************************************/
void ACEplus::processConfig(const String &filename)
{
  ACE::processConfig(filename);

  // Handle paths
  const String path=File::getPath(filename);
  char *oldPath=new char[PATH_MAX];
  getcwd(oldPath,PATH_MAX);
  chdir(path.c_str());

  // Load content sensors
  ConfigFile config(filename);
  contentSensors.setSensor(EXON,loadContentSensor("exon-model",config));
  contentSensors.setSensor(INTRON,loadContentSensor("intron-model",config));
  contentSensors.setSensor(INTERGENIC,
			   loadContentSensor("intergenic-model",config));

  // Load duration distributions
  model.intergenicLengthDistr=
    new GeometricDistribution(config.lookupOrDie("mean-intergenic-length").
			      asInt());
  model.intronLengthDistr=
    new GeometricDistribution(config.lookupOrDie("mean-intron-length").
			      asInt());
  model.exonLengthDistr=
    new EmpiricalDistribution(config.lookupOrDie("exon-length-distr"));
  model.spliceShiftDistr=
    new EmpiricalDistribution(config.lookupOrDie("splice-shift-distr"));

  // Misc initialization
  model.signalSensors=&sensors;
  model.contentSensors=&contentSensors;
  model.MAX_SPLICE_SHIFT=MAX_SPLICE_SHIFT;
  model.MIN_EXON_LEN=MIN_EXON_LEN;
  model.MIN_INTRON_LEN=MIN_INTRON_LEN;
  model.NMD_DISTANCE_PARM=NMD_DISTANCE_PARM;
  model.MAX_CRYPTIC_EXON_LEN=config.getIntOrDie("max-cryptic-exon-length");
  model.EXON_STRENGTHENING_THRESHOLD=
    config.getFloatOrDie("exon-strengthening-threshold");
  model.EXON_WEAKENING_THRESHOLD=
    config.getFloatOrDie("exon-weakening-threshold");
  model.MIN_EXON_INTRON_RATIO=
    config.getFloatOrDie("min-exon-definition-score");
  model.allowExonSkipping=allowExonSkipping;
  model.allowIntronRetention=allowIntronRetention;
  model.allowCrypticSites=allowCrypticSites;
  model.allowDeNovoSites=config.getBoolOrDie("allow-denovo-sites");
  model.allowCrypticExons=config.getBoolOrDie("allow-cryptic-exons");
  model.allowRegulatoryChanges=config.getBoolOrDie("allow-regulatory-changes");

  chdir(oldPath);
  delete [] oldPath;
}



/****************************************************************
 ACEplus::loadContentSensor()
 ****************************************************************/
ContentSensor *ACEplus::loadContentSensor(const String &label,
				       ConfigFile &config)
{
  String filename=config.lookupOrDie(label);
  return ContentSensor::load(filename);
}



/****************************************************************
 ACEplus::buildPSAs()
 ****************************************************************/
void ACEplus::buildPSAs()
{
  buildPSA(EXON);
  buildPSA(INTRON);
  buildPSA(INTERGENIC);
}



/****************************************************************
 ACEplus::buildPSA()
 ****************************************************************/
void ACEplus::buildPSA(ContentType type)
{
  PrefixSumArray &psa=contentSensors.getPSA(type);
  psa.resize(altSeqLen);
  ContentSensor *sensor=contentSensors.getSensor(type);
  psa.compute(*sensor,altSeq,altSeqStr);
}



/****************************************************************
 ACEplus::getRefLikelihood()
 ****************************************************************/
double ACEplus::getRefLikelihood(const Labeling &refLab,GffTranscript *altTrans)
{
  ProjectionChecker checker(*altTrans,*refTrans,altSeqStr,altSeq,
			    refSeqStr,refSeq,refLab,sensors);
  TranscriptSignals *signals=checker.findBrokenSpliceSites();
  GraphBuilder graphBuilder(*refTrans,*signals,model,altSeq,altSeqStr,
			    refSeq,refSeqStr,*alignment);
  LightGraph *G=graphBuilder.getGraph();
  TranscriptPaths paths(*G);
  if(paths.numPaths()!=1)
    throw String("Wrong number of reference paths: ")+paths.numPaths();
  TranscriptPath *path=paths[0];
  path->computeScore();
  double score=path->getScore();
  delete signals; delete G;
  return score;
}



/****************************************************************
 ACEplus::checkProjection()
 ****************************************************************/
void ACEplus::checkProjection(const String &outGff,bool &mapped,
			   const Labeling &refLab,
			   const Labeling &projectedLab,
			   ostream &osACE)
{
  // Convert to essex and add variants
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

  // Decompose transcript into signals
  ProjectionChecker checker(*refTrans,*altTrans,refSeqStr,refSeq,
			    altSeqStr,altSeq,projectedLab,sensors);
  String altProtein;
  if(refTrans->isCoding())
    checker.translate(*refTrans,*altTrans,refProtein,altProtein);
  TranscriptSignals *signals=checker.findBrokenSpliceSites();
  if(!signals) {
    status->prepend("unequal-numbers-of-exons");
    delete altTrans;
    return; }
  bool isCoding=altTrans->numExons()>0;
  if(isCoding) {
    int start=altTrans->getIthExon(0).getBegin();
    signals->setStartCodon(start);
  }

  // Build graph
  GraphBuilder graphBuilder(*altTrans,*signals,model,refSeq,refSeqStr,
			    altSeq,altSeqStr,*revAlignment);
  LightGraph *G=graphBuilder.getGraph();
  //cout<<*G<<endl;

  // Extract paths
  TranscriptPaths paths(*G);
  //if(graphBuilder.mapped()) { // mapped = no structure changes
  if(graphBuilder.mapped() && paths.numPaths()==1) {
    if(signals->anyWeakened()) appendBrokenSignals(signals);
    status->prepend("mapped");
    mapped=true;
    if(refTrans->isCoding()) 
      handleCoding(altTrans,checker,projectedLab);
    else
      handleNoncoding(altTrans);
  }
  else { // Enumerate alternative structures
    enumerateAlts(paths,altTransEssex,signals,altTrans,osACE,refLab,
		  projectedLab);
  }
  delete G;
  delete signals;
  delete altTrans;
}



/****************************************************************
 ACEplus::enumerateAlts()
 ****************************************************************/
void ACEplus::enumerateAlts(TranscriptPaths &paths,
			 Essex::CompositeNode *altTransEssex,
			 TranscriptSignals *signals,
			 GffTranscript *altTrans,
			 ostream &osACE,
			 const Labeling &refLab,
			 const Labeling &projectedLab)
{
  altTransEssex->deleteChild("translation");
  appendBrokenSignals(signals);
  int numPaths=paths.numPaths();
  //cout<<numPaths<<" paths"<<endl;
  if(numPaths==0) {
    status->prepend("no-transcript");
    return; }
  status->prepend("splicing-changes");
  Essex::CompositeNode *altStructNode=
    new Essex::CompositeNode("alternate-structures");
  status->append(altStructNode);

  // Compute posteriors
  //paths.computePosteriors();
  const double refLik=getRefLikelihood(refLab,altTrans);
  paths.computeLRs(refLik);

  // Process each path
  for(int i=0 ; i<numPaths ; ++i) {
    TranscriptPath *path=paths[i];
    processAltStructure(*path,altStructNode,projectedLab,i,*signals);
  }
}



void ACEplus::processAltStructure(TranscriptPath &path,
			       Essex::CompositeNode *altStructNode,
			       const Labeling &projectedLab,
			       int whichStructure,
			       TranscriptSignals &signals)
{
  // Make a GffTranscript object
  String transcriptID=
    String("ALT")+whichStructure+"_"+refTrans->getTranscriptId();
  GffTranscript *transcript=path.toTranscript(transcriptID,
					      refTrans->getGeneId(),
					      substrate,refTrans->getStrand(),
					      "ACEplus");

  // Identify reading frame
  Essex::CompositeNode *startMsg=NULL;
  if(refTrans->isCoding())
    refineStartCodon(signals.getStartCodon(),*transcript,startMsg);

  // Check for NMD
  transcript->loadSequence(altSeqStr);
  int ejcDistance;
  ProteinFate fate=nmd.predict(*transcript,altSeqStr,ejcDistance);

  // Convert to AlternativeStructure object
  AlternativeStructure alt(transcript,fate);
  alt.structureChange=path.getChange();
  Vector<ACEplus_Vertex*> vertices;
  path.getVertices(vertices);
  for(Vector<ACEplus_Vertex*>::iterator cur=vertices.begin(), end=vertices.end() ;
      cur!=end ; ++cur) {
    ACEplus_Vertex *vertex=*cur;
    if(!vertex->isAnnotated()) {
      TranscriptSignal signal(vertex->getType(),vertex->getBegin(),
			      vertex->getScore());
      signal.cryptic=true;
      signal.cutoff=vertex->getThreshold();
      signal.seq=vertex->getSeq();
      alt.crypticSignals.push_back(signal);
    }
  }

  // Check for coding change, annotate with variants
  Essex::CompositeNode *node=
    processAltStructure(alt,altStructNode,projectedLab);
  if(startMsg) node->prepend(startMsg);
}



Essex::CompositeNode *ACEplus::processAltStructure(AlternativeStructure &s,
			       Essex::CompositeNode *altStructNode,
			       const Labeling &projectedLab)
{
  Essex::CompositeNode *node=
    ACE::processAltStructure(s,altStructNode,projectedLab);
  return node;
}



void ACEplus::refineStartCodon(int start,GffTranscript &transcript,
			    Essex::CompositeNode *&msg)
{
  int newStart=StartCodonFinder::findStartCodon(transcript,
						transcript.peekExons(),
						altSeqStr,
						start,
						sensors);
  if(newStart!=start)
    if(newStart>0) {
      msg=new Essex::CompositeNode("start-codon-change");
      if(transcript.getStrand()==FORWARD_STRAND) {
	msg->append("from",start);
	msg->append("to",newStart); }
      else {
	msg->append("from",altSeqLen-start);
	msg->append("to",altSeqLen-newStart); }
    }
    else msg=new Essex::CompositeNode("start-codon-lost");
  if(newStart>=0)
    transcript.splitUTRandCDS(altSeqStr,newStart,sensors.stopCodons);
}


