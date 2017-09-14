/****************************************************************
 Mutate.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/Array2D.H"
#include "BOOM/Exceptions.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/Environment.H"
#include "BOOM/GffReader.H"
#include "Mutate.H"
using namespace std;
using namespace BOOM;

/****************************************************************
 Globals
 ****************************************************************/
static bool VERBOSE=true;
static const char *PROGRAM_NAME="Mutate";
static const char *VERSION="1.0";



/****************************************************************
 AllenMatrix::AllenMatrix()
 ****************************************************************/
AllenMatrix::AllenMatrix()
{
  matrix.resize(64,4);
}



/****************************************************************
 AllenMatrix::load()
 ****************************************************************/
void AllenMatrix::load(const String &filename)
{
  File IN(filename,"r");
  IN.getline(); // header
  while(!IN.eof()) {
    String line=IN.getline();
    Vector<String> fields;
    line.getFields(fields);
    if(fields.size()<5) continue;
    Sequence seq(fields[0],PureDnaAlphabet::global());
    int row=seq.asInt(4); // convert to a base-4 number for indexing
    float sum=0;
    for(int i=0 ; i<4 ; ++i) sum+=matrix[row][i]=fields[i+1].asFloat();
    for(int i=0 ; i<4 ; ++i) matrix[row][i]/=sum;
  }
  prepareMultinomials();
}



/****************************************************************
 AllenMatrix::sample()
 ****************************************************************/
Symbol AllenMatrix::sample(int pos,const Sequence &chrom)
{
  Sequence seq;
  if(pos<1 || pos>=chrom.getLength()-1)
    throw RootException("AllenMatrix::sample() : no context");
  chrom.getSubsequence(pos-1,3,seq);
  for(int i=0 ; i<3 ; ++i) if(seq[i]==Symbol(4)) seq[i]=3; // account for N's
  const int row=seq.asInt(4);
  return Symbol(multinomials[row].spin());
}



/****************************************************************
 AllenMatrix::prepareMultinomials()
 ****************************************************************/
void AllenMatrix::prepareMultinomials()
{
  multinomials.resize(64);
  for(int from=0 ; from<64 ; ++from) {
    RouletteWheel &wheel=multinomials[from];
    Array2D<float>::RowIn2DArray<float> row=matrix[from];
    for(int to=0 ; to<4 ; ++to) wheel.addSector(row[to]);
    wheel.doneAddingSectors();
  }
}



/****************************************************************
 SubstMatrix::prepareMultinomials()
 ****************************************************************/
void SubstMatrix::prepareMultinomials()
{
  for(int from=0 ; from<4 ; ++from) {
    RouletteWheel &wheel=multinomials[from];
    Array2D<float>::RowIn2DArray<float> row=matrix[from];
    for(int to=0 ; to<4 ; ++to) wheel.addSector(row[to]);
    wheel.doneAddingSectors();
  }
}



/****************************************************************
 SubstMatrix::sample()
 ****************************************************************/
Symbol SubstMatrix::sample(Symbol from)
{
  return Symbol(multinomials[from].spin());
}



/****************************************************************
 Mutate::Mutate()
 ****************************************************************/
Mutate::Mutate()
  : totalMutations(0), brokenSites(0), denovoDonors(0), denovoAcceptors(0)
{
  // ctor

  randomize();
  bin=Environment::lookup("ACEPLUS");
  gffTempFile=TempFilename::get();
  refTempFile=TempFilename::get();
  altTempFile=TempFilename::get();
  outTempFile=TempFilename::get();
  essexTempFile=TempFilename::get();
}



/****************************************************************
 Mutate::main()
 ****************************************************************/
int Mutate::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  parseCommandLine(cmd);
  processConfig(configFile);

  // Load substitution matrices
  //loadSubstMatrices();
  matrix.load(matrixFile);

  // Load GFF
  GffReader gffReader(refGffFile);
  Vector<GffTranscript*> *transcripts=gffReader.loadTranscripts();
  Map<String,Vector<GffTranscript*> > bySubstrate;
  for(Vector<GffTranscript*>::iterator cur=transcripts->begin(), end=
	transcripts->end() ; cur!=end ; ++cur) {
    GffTranscript *transcript=*cur;
    String substrate=transcript->getGeneId()+"_1";
    bySubstrate[substrate].push_back(transcript);
  }

  // Process FASTA file
  FastaReader reader(refFasta);
  String def, seq, id, extra;
  int genesProcessed=0;
  Regex coordRegex("(\\S+):(\\S+)-(\\S+):(\\S)");
  while(reader.nextSequence(def,seq)) {
    reader.parseDefline(def,id,extra);
    if(!bySubstrate.isDefined(id)) continue;
    Vector<GffTranscript*> &transcripts=bySubstrate[id];
    if(transcripts.size()==0) continue;
    Map<String,String> attributes;
    FastaReader::parseAttributes(extra,attributes);
    if(!attributes.isDefined("coord")) throw "no /coord on defline";
    if(!coordRegex.match(attributes["coord"]))
      throw RootException(String("can't parse /coord: "+attributes["coord"]));
    String substrate=coordRegex[1];
    int begin=coordRegex[2].asInt(), end=coordRegex[3].asInt();
    //const SubstMatrix &matrix=getMatrix(substrate,begin,end);
    Strand strand=transcripts[0]->getStrand();
    const int L=seq.length();
    if(strand==REVERSE_STRAND) {
      seq=ProteinTrans::reverseComplement(seq);
      for(Vector<GffTranscript*>::iterator cur=transcripts.begin(), end=
	    transcripts.end() ; cur!=end ; ++cur)
	(*cur)->reverseComplement(L);
    }
    Vector<SpliceSite> spliceSites;
    getSpliceSites(transcripts,spliceSites,seq);
    simulate(seq,spliceSites,matrix,transcripts);
    ++genesProcessed;
    if(genesProcessed>=MAX_GENES) break;
  }

  cerr<<"[done]"<<endl;
  return 0;
}



/****************************************************************
 Mutate::getSpliceSites()
 ****************************************************************/
void Mutate::getSpliceSites(const Vector<GffTranscript*> &transcripts,
			    Vector<SpliceSite> &into,const String &seq)
{
  // Precondition: transcripts have been rendered onto the forward strand

  Set<int> seen;
  for(Vector<GffTranscript*>::const_iterator cur=transcripts.begin(), end=
	transcripts.end() ; cur!=end ; ++cur) {
    const GffTranscript &transcript=**cur;
    Vector<Interval> introns;
    transcript.getIntrons(introns);
    for(Vector<Interval>::const_iterator cur=introns.begin(), end=
	  introns.end() ; cur!=end ; ++cur) {
      const Interval &intron=*cur;
      int donorPos=intron.getBegin();
      int acceptorPos=intron.getEnd()-2;
      if(!seen.isMember(donorPos)) {
	String consensus=seq.substring(donorPos,2);
	if(sensors.donorConsensuses.isMember(consensus))
	  into.push_back(SpliceSite(GT,donorPos));
      }
      if(!seen.isMember(acceptorPos)) {
	String consensus=seq.substring(acceptorPos,2);
	if(sensors.acceptorConsensuses.isMember(consensus))
	  into.push_back(SpliceSite(AG,acceptorPos));
      }
      seen+=donorPos; seen+=acceptorPos;
    }
  }
}



/****************************************************************
 Mutate::loadInputs()
 ****************************************************************/
void Mutate::loadInputs(const String &configFile,const String &refGffFile,
			const String &refFasta)
{
  processConfig(configFile);
  refSeqStr=loadSeq(refFasta);
  refSeq.copyFrom(refSeqStr,alphabet);
  refSeqLen=refSeqStr.length();
  refTrans=loadGff(refGffFile);
  refTrans->loadSequence(refSeqStr);
}



/****************************************************************
 Mutate::parseCommandLine()
 ****************************************************************/
void Mutate::parseCommandLine(const CommandLine &cmd)
{
  if(cmd.numArgs()!=7)
    throw String("\n\
mutate <ace.config> <ref.gff> <ref.fasta> <max-genes> <mutations-per-gene> <subst-matrices> <min-score>\n\
\n");
  configFile=cmd.arg(0);
  refGffFile=cmd.arg(1);
  refFasta=cmd.arg(2);
  MAX_GENES=cmd.arg(3).asInt();
  MUTATIONS_PER_GENE=cmd.arg(4).asInt();
  matrixFile=cmd.arg(5);
  MIN_SCORE=cmd.arg(6).asFloat();
}



/****************************************************************
 Mutate::processConfig()
 ****************************************************************/
void Mutate::processConfig(const String &filename)
{
  // Handle paths
  const String path=File::getPath(filename);
  char *oldPath=new char[PATH_MAX];
  getcwd(oldPath,PATH_MAX);
  chdir(path.c_str());

  // Load the config file
  ConfigFile config(filename);

  // Load signal sensors
  parseConsensusList("donor-consensus",config,sensors.donorConsensuses);
  parseConsensusList("acceptor-consensus",config,sensors.acceptorConsensuses);
  sensors.donorSensor=loadModel("donor-model",config);
  sensors.acceptorSensor=loadModel("acceptor-model",config);
  sensors.setConsensuses();

  // Misc initialization
  model.signalSensors=&sensors;

  // Use LLR for splice site signal sensors
  if(config.isDefined("splice-background-model")) {
    ContentSensor *bg=loadContentSensor("splice-background-model",config);
    model.contentSensors->setSpliceBackgroundModel(bg);}

  chdir(oldPath);
  delete [] oldPath;
}



/****************************************************************
 Mutate::loadContentSensor()
 ****************************************************************/
ContentSensor *Mutate::loadContentSensor(const String &label,
				       ConfigFile &config)
{
  String filename=config.lookupOrDie(label);
  return ContentSensor::load(filename);
}



/****************************************************************
 Mutate::loadSubstMatrices()
 ****************************************************************/
/*void Mutate::loadSubstMatrices()
{
  Regex regex("(\\S)->(\\S)=(\\S+)");
  File IN(matrixFile,"r");
  while(!IN.eof()) {
    String line=IN.getline();
    Vector<String> fields;
    line.getFields(fields);
    if(fields.size()<19) continue;
    SubstMatrix M;
    String substrate=M.substrate=fields[0];
    M.interval=Interval(fields[1].asInt(),fields[2].asInt());
    Array2D<float> &matrix=M.matrix;
    for(int i=3 ; i<19 ; ++i) {
      if(!regex.match(fields[i]))
	throw RootException(String("Can't parse ")+fields[i]);
      Symbol from=PureDnaAlphabet::global().lookup(regex[1][0]);
      Symbol to=PureDnaAlphabet::global().lookup(regex[2][0]);
      float P=regex[3].asFloat();
      matrix[from][to]=P;
    }
    for(int from=0 ; from<4 ; ++from) {
      float sum=0.0;
      for(int to=0 ; to<4 ; ++to) sum+=matrix[from][to];
      for(int to=0 ; to<4 ; ++to) matrix[from][to]/=sum;
    }
    matrices[substrate].push_back(M);
    matrices[substrate][matrices[substrate].size()-1].prepareMultinomials();
  }
  }*/



/****************************************************************
 Mutate::getMatrix()
 ****************************************************************/
 /*const SubstMatrix &Mutate::getMatrix(const String &substrate,int begin,
				     int end)
{
  Interval I(begin,end);
  if(!matrices.isDefined(substrate))
    throw RootException(String("no matrix found for substrate ")+substrate);
  Vector<SubstMatrix> &mats=matrices[substrate];
  int bestOverlap=-1;
  const SubstMatrix *best=NULL;
  for(Vector<SubstMatrix>::const_iterator cur=mats.begin(), end=mats.end() ;
      cur!=end ; ++cur) {
    const SubstMatrix &M=*cur;
    int overlap=M.interval.intersect(I).length();
    if(overlap>bestOverlap) {
      best=&M;
      bestOverlap=overlap;
    }
  }
  if(!best) throw "no matrix found";
  return *best;
  }*/



/****************************************************************
 Mutate::simulate()
 ****************************************************************/
void Mutate::simulate(const String &seqStr_const,
		      const Vector<SpliceSite> &spliceSites,
		      const AllenMatrix &matrix,
		      Vector<GffTranscript*> &transcripts)
{
  String seqStr=seqStr_const;
  Sequence seq(seqStr,alphabet);
  Symbol N(3), T(4);
  const int L=seqStr.length();
  for(int i=0 ; i<MUTATIONS_PER_GENE ; ++i) {
    const int pos=RandomNumber(L-2)+1;
    char oldChar=seqStr[pos];
    Symbol oldSymbol=seq[pos];
    if(oldSymbol==T) oldSymbol=N; // account for N's
    //Symbol newSymbol=matrix.sample(oldSymbol);
    Symbol newSymbol=matrix.sample(pos,seq);
    if(newSymbol==N) newSymbol=T; // account for N's
    char newChar=alphabet.lookup(newSymbol);
    seqStr[pos]=newChar;
    seq[pos]=newSymbol;
    
    // Check for broken splice sites
    bool report=false;
    if(checkBroken(seq,seqStr,spliceSites,pos)) report=true;

    // Check for de novo splice sites
    if(checkDenovo(seq,seqStr,seqStr_const,pos,spliceSites,transcripts))
      report=true;
    if(report) {
      int totalDenovo=denovoDonors+denovoAcceptors;
      cout<<totalMutations<<" mutations, "<<brokenSites<<" broken, "
	  <<totalDenovo<<" denovo ("<<denovoDonors<<" GT, "
	  <<denovoAcceptors<<" AG)"<<endl;
    }

    // Revert to the original sequence
    seqStr[pos]=oldChar;
    seq[pos]=oldSymbol;
    ++totalMutations;
  }
}



/****************************************************************
 Mutate::checkBroken_old()
 ****************************************************************/
bool Mutate::checkBroken_old(const Sequence &seq,const String &seqStr,
			 const Vector<SpliceSite> &spliceSites,
			 const int mutationPos)
{
  for(Vector<SpliceSite>::const_iterator cur=spliceSites.begin(), 
	end=spliceSites.end() ; cur!=end ; ++cur) {
    const SpliceSite &site=*cur;
    int ssPos=site.pos;
    if(mutationPos==ssPos || mutationPos==ssPos+1) {
      String consensus=seqStr.substring(ssPos,2);
      if(!isConsensus(site.type,consensus)) {
	++brokenSites;
	return true;
      }
      return false;
    }
  }
  return false;
}



/****************************************************************
 Mutate::checkBroken()
 ****************************************************************/
bool Mutate::checkBroken(const Sequence &seq,const String &seqStr,
			 const Vector<SpliceSite> &spliceSites,
			 const int mutationPos)
{
  for(Vector<SpliceSite>::const_iterator cur=spliceSites.begin(), 
	end=spliceSites.end() ; cur!=end ; ++cur) {
    const SpliceSite &site=*cur;
    SignalSensor *sensor=sensors.findSensor(site.type);
    int offset=sensor->getConsensusOffset();
    int windowLen=sensor->getContextWindowLength();
    int ssPos=site.pos;
    int windowBegin=ssPos-offset;
    int windowEnd=windowBegin+windowLen;
    if(mutationPos>=windowBegin && mutationPos<windowEnd) {
      String consensus=seqStr.substring(ssPos,2);
      if(!isConsensus(site.type,consensus)) {
	++brokenSites;
	return true;}
      else if(sensor->getLogP(seq,seqStr,windowBegin)<sensor->getCutoff()) {
	++brokenSites;
	return true;}
      return false;
    }
  }
  return false;
}



/****************************************************************
 Mutate::isConsensus()
 ****************************************************************/
bool Mutate::isConsensus(SignalType t,const String &s)
{
  switch(t) {
  case GT: return sensors.donorConsensuses.isMember(s);
  case AG: return sensors.acceptorConsensuses.isMember(s);
  default: INTERNAL_ERROR
  }
}



/****************************************************************
 Mutate::hitExistingSite()
 ****************************************************************/
bool Mutate::hitsExistingSite(const int mutationPos,
			      const Vector<SpliceSite> &spliceSites)
{
  for(Vector<SpliceSite>::const_iterator cur=spliceSites.begin(), 
	end=spliceSites.end() ; cur!=end ; ++cur) {
    const SpliceSite &site=*cur;
    int ssPos=site.pos;
    if(mutationPos==ssPos || mutationPos==ssPos+1) return true;
  }
  return false;
}



/****************************************************************
 Mutate::createsConsensus()
 ****************************************************************/
bool Mutate::createsConsensus(const String &oldSeq,const String &newSeq,
			      int pos,SignalType &t)
{
  if(pos<0 || pos+1>=oldSeq.length()) return false;
  const String oldConsensus=oldSeq.substring(pos,2);
  if(isConsensus(GT,oldConsensus) || isConsensus(AG,oldConsensus))
    return false;
  const String newConsensus=newSeq.substring(pos,2);
  bool createsGT=false, createsAG=false;
  if(isConsensus(GT,newConsensus)) createsGT=true;
  if(isConsensus(AG,newConsensus)) createsAG=true;
  if(createsGT && createsAG) {
    if(Random0to1()<0.5) t=GT;
    else t=AG; }
  else if(createsGT) t=GT;
  else if(createsAG) t=AG;
  return createsGT || createsAG;
}



/****************************************************************
 Mutate::incrementDenovo()
 ****************************************************************/
void Mutate::incrementDenovo(SignalType t)
{
  if(t==GT) ++denovoDonors;
  else if(t==AG) ++denovoAcceptors;
  else INTERNAL_ERROR
}



/****************************************************************
 Mutate::createsSite()
 ****************************************************************/
bool Mutate::createsSite(const Sequence &seq,const String &seqStr,
			 int consensusPos,SignalType t)
{
  SignalSensor *sensor=sensors.findSensor(t);
  const float threshold=sensor->getCutoff();
  const int offset=sensor->getConsensusOffset();
  const int windowLen=sensor->getContextWindowLength();
  const int begin=consensusPos-offset;
  const int end=begin+windowLen;
  if(begin<0 || end>seqStr.length()) return false;
  const float logP=sensor->getLogP(seq,seqStr,begin);
  return logP>=threshold;
}



/****************************************************************
 Mutate::checkDenovo()
 ****************************************************************/
bool Mutate::checkDenovo(const Sequence &seq,const String &newSeqStr,
			 const String &oldSeqStr,const int mutationPos,
			 const Vector<SpliceSite> &spliceSites,
			 Vector<GffTranscript*> &transcripts)
{
  if(hitsExistingSite(mutationPos,spliceSites)) return false;
  SignalType t;
  if(createsConsensus(oldSeqStr,newSeqStr,mutationPos,t)) {
    if(createsSite(seq,newSeqStr,mutationPos,t)) {
      incrementDenovo(t);
      if(MIN_SCORE<=0.0 ||
	 scoreDenovo(seq,newSeqStr,oldSeqStr,transcripts)>=MIN_SCORE)
	return true; }
  }
  if(createsConsensus(oldSeqStr,newSeqStr,mutationPos-1,t)) {
    if(createsSite(seq,newSeqStr,mutationPos-1,t)) {
      incrementDenovo(t);
      if(MIN_SCORE<=0.0 ||
	 scoreDenovo(seq,newSeqStr,oldSeqStr,transcripts)>MIN_SCORE)
	return true; }
  }
  return false;
}



/****************************************************************
 Mutate::loadScore()
 ****************************************************************/
float Mutate::loadScore(GffTranscript &refTrans,
			Vector<GffTranscript*> &transcripts)
{
  float bestScore=NEGATIVE_INFINITY;
  for(Vector<GffTranscript*>::iterator cur=transcripts.begin(), 
	end=transcripts.end() ; cur!=end ; ++cur) {
    GffTranscript *transcript=*cur;
    if(transcript->identical(refTrans)) continue;
    if(!transcript->hasValidScore()) continue;
    float score=transcript->getScore();
    //cout<<"SCORE="<<score<<endl;
    if(score>bestScore) bestScore=score;
  }
  return bestScore;
}



/****************************************************************
 Mutate::scoreDenovo()
 ****************************************************************/
float Mutate::scoreDenovo(const Sequence &seq,const String &newSeqStr,
			  const String &oldSeqStr,
			  Vector<GffTranscript*> &transcripts)
{
  float bestScore=NEGATIVE_INFINITY;
  for(Vector<GffTranscript*>::iterator cur=transcripts.begin(), end=
	transcripts.end() ; cur!=end ; ++cur) {
    GffTranscript *transcript=*cur;

    // Write GFF file
    ofstream os(gffTempFile.c_str());
    transcript->toGff(os);
    os.close();

    // Write fasta files
    const int L=newSeqStr.length();
    String cigar=String("")+L+"M";
    String altDefline=String(">chr /cigar=")+cigar;
    fastaWriter.writeFasta(">chr",oldSeqStr,refTempFile);
    fastaWriter.writeFasta(altDefline,newSeqStr,altTempFile);

    // Run ACE+
    String cmd=bin+"/aceplus "+configFile+" "+gffTempFile+" "+
      refTempFile+" "+altTempFile+" "+outTempFile+" "+essexTempFile+
      " 1> /dev/null 2> /dev/null";
    system(cmd.c_str());
    String cmd2=bin+"/essex-to-gff-AS2.pl "+essexTempFile+" "+outTempFile+" 1";
    system(cmd2.c_str());

    // Parse output to get score
    GffReader reader(outTempFile);
    Vector<GffTranscript*> *predictions=reader.loadTranscripts();
    //system((String("cat ")+outTempFile).c_str());
    float score=loadScore(*transcript,*predictions);
    //cout<<"score\t"<<score<<endl;
    for(Vector<GffTranscript*>::iterator cur=predictions->begin(), 
	  end=predictions->end() ; cur!=end ; ++cur) delete *cur;
    delete predictions;
    if(score>bestScore) bestScore=score;

  /*
aceplus <ace.config> <ref.gff> <ref.fasta> <alt.fasta> <out.gff> <out.essex>
     -c = sequence has been reversed, but cigar string has not
  alt.fasta must have a cigar string: >ID ... /cigar=1045M3I10M7D4023M ...
   */
  }
  return bestScore;
}




