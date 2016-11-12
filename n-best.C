/****************************************************************
 n-best.C : find N best predictions (isoforms)
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "BOOM/Constants.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Array1D.H"
#include "BOOM/Vector.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/FixedSizePriorityQueue.H"
#include "BOOM/Time.H"
#include "BOOM/Stack.H"
#include "BOOM/WigBinary.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Sequence.H"
#include "BOOM/String.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/Set.H"
#include "LightGraph.H"
#include "TrellisLink.H"
#include "NMD.H"
#include "CIA.H"
using namespace std;
using namespace BOOM;

#ifdef EXPLICIT_GRAPHS
#error Please edit the file genezilla.H, comment out the definition of EXPLICIT_GRAPHS, issue a "make clean", and recompile this project
#endif

static const char *PROGRAM_NAME="CIA";
static const char *VERSION="0.1";

Alphabet alphabet;
int frame; // ### CAUTION: this is required by older code; to be removed

class Application {
public:
  Application();
  void main(int argc,char *argv[]);
  void buildTrellis(LightGraph &,int N,BOOM::Vector<TrellisLink> &);
protected:
  String graphFile;
  Set<String> stopCodons;
  DnaAlphabet alphabet;
  int autoNmode;
  Array1D< Array1D<TrellisLink> > links;
  Vector<LightVertex*> ATGs, annotatedATGs;
  int queueCapacity;
  TrellisLinkComparator cmp;
  bool wantIntrons;
  WigBinary *wig;
  Vector<WigInterval> nonzeroRegions;
  Array1D<bool> edgeEmitted;
  int numSupportedIntrons, supportedIntronsEmitted;
  int numSupportedExons, supportedExonsEmitted;
  int transcriptNum;
  String substrate;
  Sequence substrateSeq;
  NMD nmd;
  SignalSensor *wideATG;
  Vector<int> scannedATGs;
  int posmod(int x);
  int countSupportedIntrons(LightGraph &);
  int countSupportedExons(LightGraph &);
  void traceback(TrellisLink *endLink,Vector<TrellisLink*> &path);
  void emitHeader();
  bool generateTranslations(Vector<TrellisLink*> &,LightGraph &,int &parseNum);
  bool emit(int startPos,Vector<TrellisLink*> &,LightGraph &,int parseNum);
  void emitFooter(LightGraph &G,BOOM::Vector<TrellisLink> &parseList);
  void getExonEdges(const Vector<TrellisLink*> &links,
		    Vector<LightEdge*> &into);
  int mapToTranscriptCoords(int genomicCoord,Vector<LightEdge*> &exons);
  void getTranscriptSeq(const Vector<LightEdge*> &exons,const Sequence
			&chrSeq,const String &chrStr,Sequence &intoSeq,
			String &intoStr);
  void makeTranscript(Vector<TrellisLink*> &path,GffTranscript &);
  int findStopCodon(const String &splicedTranscript,int startPos);
  void setCDS(GffTranscript &,int cdsBegin,int cdsEnd);
  float getGCcontent(const BOOM::String &seq);
  void scanATGs();
};



int main(int argc,char *argv[]) {
  try {	Application app;  app.main(argc,argv); }
  catch(const char *p) { cerr << p << endl; return -1; }
  catch(const string &msg) { cerr << msg.c_str() << endl; return -1; }
  catch(const exception &e) {
    cerr << "STL exception caught in main:\n" << e.what() << endl;
    return -1;
  }
  catch(...) {
    cerr << "Unknown exception caught in main" << endl;
    return -1;
  }
  return 0;
}



Application::Application()
  : wantIntrons(false), queueCapacity(10), wig(NULL), transcriptNum(0)
{
  stopCodons+="TGA";
  stopCodons+="TAG";
  stopCodons+="TAA";
}



void Application::main(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"q:");
  if(cmd.numArgs()!=4)
    throw BOOM::String("\n\
n-best [options] <model.iso> <in.graph> <N> substrate.fasta\n\
    -q queue-capacity\n\
");
  const String isochoreFilename=cmd.arg(0);
  graphFile=cmd.arg(1);
  const int N=cmd.arg(2).asInt();
  const String substrateFile=cmd.arg(3);
  String def;
  FastaReader::load(substrateFile,def,substrate);
  substrateSeq=Sequence(substrate,alphabet);
  cout<<"Processing graph file: "<<graphFile<<endl;

  // Some global initialization
  alphabet=DnaAlphabet::global();
  if(cmd.option('q')) queueCapacity=cmd.optParm('q').asInt();
  wantIntrons=cmd.option('i');
  autoNmode=0;
  GarbageIgnorer gc;
  IsochoreTable isochores(gc);
  isochores.load(isochoreFilename);
  const float gcContent=getGCcontent(substrate);
  Isochore *isochore=isochores.getIsochore(gcContent);
  wideATG=isochore->wideATG;
  if(!wideATG)
    throw String("wide-start-codon-model is not defined in ")+isochoreFilename;
  scanATGs();

  // Load graph
  LightGraph *G=new LightGraph(graphFile);
  const int numEdges=G->getNumEdges();
  edgeEmitted.resize(G->getLargestEdgeID()+1);//numEdges);
  edgeEmitted.setAllTo(false);
  G->getATGs(ATGs); G->getAnnotatedATGs(annotatedATGs);

  // Build trellis
  BOOM::Vector<TrellisLink> termini;
  buildTrellis(*G,N,termini);

  // Construct paths & emit GFF
  emitHeader();
  int parseNum=1;
  const int numTermini=termini.size();
  for(int i=0 ; i<numTermini ; ++i) {
    Vector<TrellisLink*> path;

    // Generate next splice pattern
    traceback(&termini[i],path);

    // Generate zero or more translation frames & emit them as predictions
    if(!generateTranslations(path,*G,parseNum)) break;
  }
  emitFooter(*G,termini);
}



int Application::countSupportedIntrons(LightGraph &G)
{
  int n=0;
  int numEdges=G.getNumEdges();
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *e=G.getEdge(i);
    if(dropStrand(e->getEdgeType())==INTRON && e->isSupported()) ++n;
  }
  return n;
}



int Application::countSupportedExons(LightGraph &G)
{
  int n=0;
  int numEdges=G.getNumEdges();
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *e=G.getEdge(i);
    if(isCoding(dropStrand(e->getEdgeType())) && e->isSupported()) ++n;
  }
  return n;
}



int Application::posmod(int x) 
{
  int f=x%3;
  if(f>=0) return f;
  return f+3;
}



void Application::buildTrellis(LightGraph &G,int N,
			       BOOM::Vector<TrellisLink> &termini)
{
  const int numVertices=G.getNumVertices();
  links.resize(numVertices);
  if(numVertices==0) return;
  for(int i=0 ; i<numVertices ; ++i) {
    FixedSizePriorityQueue<TrellisLink> *Q[3];
    for(int j=0 ; j<3 ; ++j)
      Q[j]=new FixedSizePriorityQueue<TrellisLink>(queueCapacity,cmp);
    LightVertex *signal=G.getVertex(i);
    if(!signal) continue;
    SignalType signalType=dropStrand(signal->getSignalType());
    Strand signalStrand=signal->getStrand();
    int defaultSignalPhase=signalStrand==FORWARD_STRAND ? 0 : 2;
    Vector<LightEdge*> &edges=signal->getEdgesIn();
    int numEdges=edges.size();
    if(numEdges==0) { //left-terminus:
      TrellisLink link(NULL,NULL);
      Q[0]->insert(link);
    }
    for(Vector<LightEdge*>::iterator cur=edges.begin(), end=edges.end() ; 
	cur!=end ; ++cur) {
      LightEdge *currentEdge=*cur;
      LightVertex *pred=currentEdge->getLeft();
      Array1D<TrellisLink> &predLinks=links[pred->getID()];
      const int numPredLinks=predLinks.size();
      if(currentEdge->isExon()) { //exon (coding or UTR)
	for(int j=0 ; j<numPredLinks ; ++j) {
	  TrellisLink &predLink=predLinks[j];
	  const int phase=predLink.getPhase();
	  const int signalPhase=0;//currentEdge->propagateForward(phase);
	  if(signalStrand==FORWARD_STRAND) {
	    if((signalType==TAG || signalType==ATG) && signalPhase!=0)
	      continue; }
	  else // REVERSE_STRAND
	    if((signalType==TAG || signalType==ATG) && signalPhase!=2)
	      continue;
	  double score=currentEdge->getScore(phase)+
	    predLink.getScore();
	  if(isFinite(score)) {
	    TrellisLink link(&predLink,currentEdge,signalPhase,score);
	    Q[signalPhase]->insert(link);
	  }
	}
      }
      else if(currentEdge->isIntron()) { //intron
	for(int j=0 ; j<numPredLinks ; ++j) {
	  TrellisLink &predLink=predLinks[j];
	  const int phase=predLink.getPhase();
	  const double score=currentEdge->getScore(phase)+
	    predLink.getScore();
	  if(isFinite(score)) {
	    TrellisLink link(&predLink,currentEdge,phase,score);
	    Q[phase]->insert(link);
	  }	
	}
      }
      else { //intergenic
	for(int j=0 ; j<numPredLinks ; ++j) {
	  TrellisLink &predLink=predLinks[j];
	  const int phase=defaultSignalPhase;
	  const double score=currentEdge->getScore(0)+predLink.getScore();
	  if(isFinite(score)) {
	    TrellisLink link(&predLink,currentEdge,phase,score);
	    Q[phase]->insert(link);
	  }
	}
      }
    }

    // Copy selected links into signal's link set
    Array1D<TrellisLink> &linkSet=links[i];
    int size=0;
    for(int f=0 ; f<3 ; ++f) size+=Q[f]->getNumElements();
    linkSet.resize(size);
    int j=0;
    for(int f=0 ; f<3 ; ++f)
      for(FixedSizePriorityQueue<TrellisLink>::iterator cur=Q[f]->begin(), 
	    end=Q[f]->end() ; cur!=end ; ++cur)
	linkSet[j++]=*cur;
    for(int j=0 ; j<3 ; ++j) delete Q[j];
  }
  
  // Find all of the right-terminal links
  for(int rCur=numVertices-1; rCur>=0; --rCur) {
    LightVertex *rt=G.getVertex(rCur);
    if(!rt) continue;
    if(!rt->getEdgesOut().isEmpty()) break;
    Array1D<TrellisLink> &rtLinks=links[rCur];
    const int numRtLinks=rtLinks.size();
    for(int j=0 ; j<numRtLinks ; ++j)
      termini.push_back(rtLinks[j]);
  }
  
  // Return links representing the N best parses
  VectorSorter<TrellisLink> sorter(termini,cmp);
  sorter.sortDescendInPlace();
  if(N<termini.size()) termini.resize(N);
}



void Application::traceback(TrellisLink *endLink,Vector<TrellisLink*> &path) 
{
  Stack<TrellisLink*> pStack;
  TrellisLink *currentLink=endLink;
  while(currentLink) {
    pStack.push(currentLink);
    currentLink=currentLink->getPred();
  }
  while(!pStack.isEmpty()) {
    TrellisLink *pl=pStack.pop();
    path.push_back(pl);
  }
}



void Application::emitHeader()
{
  cout<<"##gff-version 2"<<endl;
  cout<<"##source-version Decoder 1.0"<<endl;
  cout<<"##date "<<getDateAndTime()<<endl;
  cout<<"##Type DNA"<<endl;
}



void Application::emitFooter(LightGraph &G,
			     BOOM::Vector<TrellisLink> &parseList) 
{
  cout<<"\n#======================================================================"<<endl;
  cout<<"\n# sequence length: "<<G.getSubstrateLength()<<endl;
}



void Application::getExonEdges(const Vector<TrellisLink*> &path,
			       Vector<LightEdge*> &into)
{
  for(Vector<TrellisLink*>::const_iterator cur=path.begin(), end=path.end() ;
      cur!=end ; ++cur) {
    TrellisLink *pl=*cur;
    LightEdge *e=pl->getEdge();
    if(e && (isCoding(e->getEdgeType()) || isUTR(e->getEdgeType())))
       into.push_back(e);
  }
}



int Application::mapToTranscriptCoords(int genomicCoord,
				       Vector<LightEdge*> &exons)
{
  int sum=0;
  for(Vector<LightEdge*>::iterator cur=exons.begin(), end=exons.end() ;
      cur!=end ; ++cur) {
    LightEdge *edge=*cur;
    if(edge->getStrand()!=FORWARD_STRAND) INTERNAL_ERROR;
    const int begin=edge->getBegin(), end=edge->getEnd(),
      exonLen=edge->getLength();
    if(genomicCoord<begin) return -1;
    if(genomicCoord<end) return sum+genomicCoord-begin;
    sum+=exonLen;
  }
  return -1;
}



void Application::getTranscriptSeq(const Vector<LightEdge*> &exons,
				   const Sequence &chrSeq,
				   const String &chrStr,
				   Sequence &intoSeq,
				   String &intoStr)
{
  intoSeq.clear();
  intoStr.clear();
  for(Vector<LightEdge*>::const_iterator cur=exons.begin(), end=exons.end() ;
      cur!=end ; ++cur) {
    LightEdge *edge=*cur;
    if(edge->getStrand()!=FORWARD_STRAND) INTERNAL_ERROR;
    const int begin=edge->getBegin(), end=edge->getEnd(),
      exonLen=edge->getLength();
    Sequence exonSeq;
    chrSeq.getSubsequence(begin,exonLen,exonSeq);
    intoSeq.append(exonSeq);
    intoStr+=chrStr.substring(begin,exonLen);
  }
}



bool Application::generateTranslations(Vector<TrellisLink*> &path,
				       LightGraph &G,int &parseNum)
{
  bool emitted=false;
  Vector<LightEdge*> exons;
  getExonEdges(path,exons);

  // Consider each ATG vertex (annotated ATG and any upstream ATG created
  // via a SNP or indel)
  for(Vector<LightVertex*>::iterator cur=ATGs.begin(), end=ATGs.end() ;
      cur!=end ; ++cur) {
    LightVertex *ATG=*cur;
    int startPos=mapToTranscriptCoords(ATG->getBegin(),exons);
    if(startPos<0) continue;

    // Emit prediction
    if(emit(startPos,path,G,parseNum)) { emitted=true; ++parseNum; }
    break; // ### for now, we're emitting only one translation
  }
  if(!emitted) // annotated ATG didn't survive, so scan for others
    for(Vector<int>::iterator cur=scannedATGs.begin(), end=scannedATGs.end() ;
      cur!=end ; ++cur){
      int genomicPos=*cur;
      int startPos=mapToTranscriptCoords(genomicPos,exons);
      if(startPos<0) continue;

      // Emit prediction
      if(emit(startPos,path,G,parseNum)) { emitted=true; ++parseNum; }
      break; // ### for now, we're emitting only one translation
  }
  return emitted;
}



void Application::makeTranscript(Vector<TrellisLink*> &path,
				 GffTranscript &transcript)
{
  transcript.setStrand(FORWARD_STRAND);
  for(Vector<TrellisLink*>::iterator cur=path.begin(), end=path.end() ;
      cur!=end ; ++cur) {
    TrellisLink *pl=*cur;
    LightEdge *e=pl->getEdge();
    if(!e) continue;
    if(e->getStrand()!=FORWARD_STRAND) INTERNAL_ERROR;
    const ContentType edgeType=dropStrand(e->getEdgeType());
    if(!::isCoding(edgeType) && !::isUTR(edgeType)) continue;
    const int startPos=e->getBegin();
    const int endPos=e->getEnd();
    const int phase=0; // ###
    int length=endPos-startPos+1;
    float score=exp(e->getScore(phase)/length)/.25;
    GffExon *exon=new GffExon(ET_UTR,startPos,endPos,transcript,true,score,
			      true,phase);
    transcript.addExon(exon);
  }
}



int Application::findStopCodon(const String &splicedTranscript,int startPos)
{
  const int L=splicedTranscript.length();
  for(int pos=startPos+3 ; pos+2<L ; pos+=3) {
    String codon=splicedTranscript.substring(pos,3);
    if(stopCodons.isMember(codon)) return pos;
  }
  return -1;
}



void Application::setCDS(GffTranscript &transcript,int cdsBegin,
			 int cdsEnd)
{
  const int genomicStart=transcript.mapToGenomicCoords(cdsBegin);
  const int genomicEnd=transcript.mapToGenomicCoords(cdsEnd);
  if(graphFile=="gi_170932551_ref_NM_032790.3_.graph") {
    cout<<"genomicStart="<<genomicStart<<" genomicEnd="<<genomicEnd<<endl;
    transcript.toGff(cout); cout<<endl;
  }
  GffTranscript newTrans;
  newTrans.setGeneId(transcript.getGeneId());
  newTrans.setTranscriptId(transcript.getTranscriptId());
  newTrans.setStrand(FORWARD_STRAND);
  int numExons=transcript.numExons();

  // Make 5' UTR exons
  int i;  bool seenEnd=false;
  for(i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript.getIthExon(i);
    if(genomicStart<exon.getBegin()) INTERNAL_ERROR;
    if(genomicStart<exon.getEnd()) {
      if(genomicStart>exon.getBegin()) {
	GffExon *newUTR=new GffExon(exon,newTrans);
	newUTR->setEnd(genomicStart);
	newUTR->changeExonType(ET_UTR);
	newTrans.addUTR(newUTR); }
      GffExon *newExon=new GffExon(exon,newTrans);
      newExon->setBegin(genomicStart);
      newExon->changeExonType(ET_EXON);
      newTrans.addExon(newExon);
      if(genomicEnd<newExon->getEnd()) {
	newExon->setEnd(genomicEnd);
	seenEnd=true;
	--i; }
      break;
    }
    else { 
      GffExon *newUTR=new GffExon(exon,newTrans);
      newUTR->changeExonType(ET_UTR);
      newTrans.addUTR(newUTR);
    }
  }
  if(i>=numExons) INTERNAL_ERROR;

  // Make CDS exons
  for(++i ; i<numExons ; ++i) {
    GffExon &exon=transcript.getIthExon(i);
    if(genomicEnd<exon.getBegin() && !seenEnd) INTERNAL_ERROR; // ###
    if(genomicEnd<exon.getEnd()) {
      if(genomicEnd>exon.getBegin() && !seenEnd) {
	GffExon *newExon=new GffExon(exon,newTrans);
	newExon->setEnd(genomicEnd);
	newExon->changeExonType(ET_EXON);
	newTrans.addExon(newExon);
      }
      GffExon *newExon=new GffExon(exon,newTrans);
      newExon->setBegin(genomicEnd);
      newExon->changeExonType(ET_UTR);
      newTrans.addUTR(newExon);
      break;
    }
    else {
      GffExon *newExon=new GffExon(exon,newTrans);
      newExon->changeExonType(ET_EXON);
      newTrans.addExon(newExon);
    }
  }

  // Make 3' UTR exons
  for(++i ; i<numExons ; ++i) {
    GffExon &exon=transcript.getIthExon(i);
    GffExon *newExon=new GffExon(exon,newTrans);
    newExon->changeExonType(ET_UTR);
    newTrans.addUTR(newExon);
  }

  transcript=newTrans;
  transcript.setExonTypes();
  transcript.setUTRtypes();
  transcript.computePhases();
}



bool Application::emit(int startPos,Vector<TrellisLink*> &path,LightGraph &G,
		       int parseNum)
{
  // First, make a basic transcript
  TrellisLink *end=path[path.size()-1];
  const float transcriptScore=exp(end->getScore()/G.getSubstrateLength())/.25;
  GffTranscript transcript(String(transcriptNum),G.getSubstrate(),
			   FORWARD_STRAND,"CIA");
  makeTranscript(path,transcript);

  // Find end of CDS (stop codon)
  transcript.loadSequence(substrate);
  String splicedTranscript=transcript.getFullSequence();
  const int splicedLen=splicedTranscript.length();
  if(startPos+3>splicedLen) return false;
  const int stopPos=findStopCodon(splicedTranscript,startPos);
  const int endOfCDS=stopPos>0 ? stopPos+3 : splicedLen;
  
  // Now edit the transcript structure by putting in the start & stop codons
  setCDS(transcript,startPos,endOfCDS);

  // Generate output
  cout<<"\n#======================================================================"<<endl;
  cout<<"\n# transcript structure #"<<parseNum<<":"<<endl;
  transcript.toGff(cout);

  String flags;
  if(substrate.length()>0) {
    const NMD_TYPE ptcType=nmd.predict(transcript,substrate);
    switch(ptcType) {
    case NMD_NONE: break;
    case NMD_NMD: flags="PTC=NMD;"; break;
    case NMD_TRUNCATION: flags="PTC=truncation;"; break;
    case NMD_NO_STOP: flags="PTC=no_stop_codon;"; break;
    case NMD_NO_START: flags="PTC=no_start_codon;"; break;
    }
  }
  const String contigID=G.getSubstrate();
  cout<<contigID << "\t" << "CIA" << "\t" << "transcript" << 
    "\t" << transcript.getBegin()+1 << "\t" << transcript.getEnd() << "\t" 
      << transcriptScore << "\t" << transcript.getStrand() << "\t" << "." 
      << "\ttranscript_id=" << transcriptNum << "; " << flags <<endl;

  return true;
}




float Application::getGCcontent(const BOOM::String &seq)
{
  int n=seq.length(), ATCG=0, GC=0;
  const char *p=seq.c_str();
  for(int i=0 ; i<n ; ++i)
    switch(*p++)
      {
      case 'G':
      case 'C':
	++GC;
	// no break...fall through:
      case 'A':
      case 'T':
	++ATCG;
      }
  return GC/float(ATCG);
}



 void Application::scanATGs()
 {
   const int L=substrate.length();
   const int windowLen=wideATG->getContextWindowLength();
   const int consensusOffset=wideATG->getConsensusOffset();
   const float threshold=wideATG->getCutoff();
   const int lastPos=L-windowLen;
   for(int i=0 ; i<=lastPos ; ++i) {
     const int consensusPos=i+consensusOffset;
     if(!wideATG->consensusOccursAt(substrate,consensusPos)) continue;
     const float logP=wideATG->getLogP(substrateSeq,substrate,i);
     if(logP>=threshold) scannedATGs.push_back(i);
   }
 }



