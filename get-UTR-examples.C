/****************************************************************
 get-UTR-examples.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/GffReader.H"
#include "BOOM/Regex.H"
#include "BOOM/Vector.H"
#include "BOOM/Interval.H"
#include "BOOM/Array2D.H"
#include "SignalType.H"
#include "ContentType.H"
using namespace std;
using namespace BOOM;

class Application {
  FastaWriter writer;
  Regex chromFilename;
  Vector<FloatInterval> isochores;
  Vector<ostream*> outfiles5, outfiles3;
  void closeFiles(Vector<ostream*> &);
  float getGC(const String &);
  ostream *getIsochoreFile(float gc,Vector<ostream*> &);
  void getTransProbs(Map<String,Vector<GffGene> > &,const String &outfile);
  void updateTransCounts(const GffTranscript &,Array2D<int> &counts);
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
  : chromFilename("(chr\\S+)\.fa")
{
  // ctor

  isochores.push_back(FloatInterval(0.0,43.0));
  isochores.push_back(FloatInterval(43.0,51.0));
  isochores.push_back(FloatInterval(51.0,57.0));
  isochores.push_back(FloatInterval(57.0,100.1));
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("get-UTR-examples <*.gff> <genome-dir> <out-filestem>");
  const String gffFile=cmd.arg(0);
  const String genomeDir=cmd.arg(1);
  const String filestem=cmd.arg(2);

  // Open output files
  for(Vector<FloatInterval>::iterator cur=isochores.begin(), end=
	isochores.end() ; cur!=end ; ++cur) {
    FloatInterval I=*cur;
    String filename5=filestem+".5prime.iso"+int(I.getBegin())+"-"+
      int(I.getEnd())+".fasta";
    outfiles5.push_back(new ofstream(filename5.c_str()));
    String filename3=filestem+".3prime.iso"+int(I.getBegin())+"-"+
      int(I.getEnd())+".fasta";
    outfiles3.push_back(new ofstream(filename3.c_str()));
  }

  // Load GFF
  Map<String,Vector<GffGene> > &genesByChrom=*GffReader::genesByChrom(gffFile);

  // Estimate transition probabilities
  const String transFile=filestem+".trans";
  getTransProbs(genesByChrom,transFile);

  // Process each chromosome separately
  Vector<String> chromFiles;
  File::getFileList(genomeDir,chromFiles);
  for(Vector<String>::iterator cur=chromFiles.begin(), end=chromFiles.end() ;
      cur!=end ; ++cur) {
    const String chromFile=*cur;
    if(!chromFilename.search(chromFile)) continue;
    const String chr=chromFilename[1];
    if(!genesByChrom.isDefined(chr)) continue;
    Vector<GffGene> &genes=genesByChrom[chr];
    String def, seq;
    FastaReader::load(genomeDir+"/"+chromFile,def,seq);
    for(Vector<GffGene>::iterator cur=genes.begin(), end=genes.end() ;
	cur!=end ; ++cur) {
      GffGene &gene=*cur;
      GffTranscript &transcript=*gene.longestTranscript();
      transcript.loadSequence(seq);
      const float gc=getGC(transcript.getSequence());
      ostream *os5=getIsochoreFile(gc,outfiles5);
      ostream *os3=getIsochoreFile(gc,outfiles3);
      transcript.setUTRtypes();
      String id=transcript.getTranscriptId();
      int exonNum=0;
      for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	    transcript.getUTRend() ; cur!=end ; ++cur) {
	GffExon *exon=*cur;
	ostream *os=exon->isUTR5() ? os5 : os3;
	const int prime=exon->isUTR5() ? 5 : 3;
	String def=String(">")+id+"_UTR"+prime+"_"+exonNum
	  +" /chr="+chr+" /begin="+exon->getBegin()+" /end="+exon->getEnd();
	writer.addToFasta(def,exon->getSequence(),*os);
	++exonNum;
      }
    }
  }

  // Clean up
  closeFiles(outfiles5);
  closeFiles(outfiles3);
  return 0;
}




float Application::getGC(const String &S)
{
  const int A=S.count('A'), C=S.count('C'), G=S.count('G'), T=S.count('T');
  const int L=A+G+C+T;
  return 100*(G+C)/float(L);
}



ostream *Application::getIsochoreFile(float gc,Vector<ostream*> &outfiles)
{
  const int N=isochores.size();
  for(int i=0 ; i<N ; ++i)
    if(isochores[i].contains(gc)) return outfiles[i];
  INTERNAL_ERROR;
}



void Application::closeFiles(Vector<ostream*> &files)
{
  for(Vector<ostream*>::iterator cur=files.begin(), end=files.end() ;
      cur!=end ; ++cur)
    delete *cur;
}



void Application::getTransProbs(Map<String,Vector<GffGene> > &genesByChrom,
				const String &outfile)
{
  Array2D<int> counts(numSignalTypes(),numSignalTypes()); counts.setAllTo(0);
  ofstream os(outfile.c_str());
  Set<String> chroms;
  genesByChrom.getKeys(chroms);
  for(Set<String>::iterator cur=chroms.begin(), end=chroms.end() ; cur!=end ;
      ++cur) {
    const String &chr=*cur;
    Vector<GffGene> &genes=genesByChrom[chr];
    for(Vector<GffGene>::iterator cur=genes.begin(), end=genes.end() ;
	cur!=end ; ++cur) {
      GffGene &gene=*cur;
      GffTranscript &transcript=*gene.longestTranscript();
      updateTransCounts(transcript,counts);
    }
  }

  // Normalize counts into probabilities and emit
  const int N=numSignalTypes();
  for(SignalType t=0 ; t<N ; ++t) {
    if(getStrand(t)==REVERSE_STRAND) continue;
    int sum=0;
    for(SignalType u=0 ; u<N ; ++u) {
      if(getStrand(u)==REVERSE_STRAND) continue;
      const int count=counts[t][u];
      sum+=count;
    }
    if(sum==0) continue;
    for(SignalType u=0 ; u<N ; ++u) {
      if(getStrand(u)==REVERSE_STRAND) continue;
      const int count=counts[t][u];
      if(count==0) continue;
      float P=float(count)/sum;
      os<<t<<" -> "<<u<<" : "<<P<<endl;
    }
  }
}



void Application::updateTransCounts(const GffTranscript &transcript,
				    Array2D<int> &counts)
{
  const int numExons=transcript.numExons()+transcript.numUTR();
  const int numIntrons=numExons-1;
  if(numExons==0) return;
  ++counts[LEFT_TERMINUS][TSS];
  ++counts[TES][RIGHT_TERMINUS];
  if(numIntrons>0) {
    ++counts[TSS][UTR5GT];
    counts[UTR5GT][UTR5AG]+=numIntrons;
    counts[UTR5AG][UTR5GT]+=numIntrons;
    ++counts[UTR5AG][TES];
  }

  /*
  // Make a vector of all exons, including UTRs, in transcription order
  Vector<GffExon*> allExons;
  Vector<GffExon*>::iterator ucur=transcript.getUTR(), 
    uend=transcript.getUTRend(), cur=transcript.getExons(),
    end=transcript.getExonsEnd();
  for(; ucur!=uend ; ++ucur) {
    GffExon *utr=*ucur;
    if(utr->isUTR3()) break;
    allExons.push_back(utr);
  }
  for(; cur!=end ; ++cur) allExons.push_back(*cur);
  for(; ucur!=uend ; ++ucur) allExons.push_back(*ucur);

  // Iterate through all exons and make a list of signals
  Vector<SignalType> signals;
  signals.push_back(LEFT_TERMINUS);
  const int numExons=allExons.size();
  for(int i=0 ; i<numExons ; ++i) {
    GffExon *exon=allExons[i];
    if(i==0) signals.push_back(TSS);
    else signals.push_back(UTR5AG);
    if(i==numExons-1) signals.push_back(TSS);
    else signals.push_back(UTR5GT);
  }
  signals.push_back(RIGHT_TERMINUS);

  // Iterate through all signals, counting transitions
  const int N=signals.size();
  for(int i=0 ; i<N-1 ; ++i) {
    const SignalType thisSignal=signals[i], nextSignal=signals[i+1];
    ++counts[thisSignal][nextSignal];
  }
  */
}




