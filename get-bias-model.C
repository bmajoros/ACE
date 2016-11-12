/****************************************************************
 get-bias-model.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "BOOM/Constants.H"
#include "BOOM/CommandLine.H"
#include "BOOM/WigBinary.H"
#include "BOOM/Array1D.H"
#include "BOOM/GffReader.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Regex.H"
#include "BOOM/StringMap.H"
using namespace std;
using namespace BOOM;


const int HASH_SIZE=9901;


class Application {
public:
  Application();
  void main(int argc,char *argv[]);
protected:
  int wordLen;
  Regex filenameRegex;
  StringMap<float> counts, biasedCounts;
  void processChunk(const String &dir,const String &stem);
  int getLongest(Vector<String> &seqs);
  void getTranscriptSeqs(const String &dir,const String &stem,
			 Vector<String> &seqs,
			 Vector<GffTranscript*> &transcripts);
  void updateStats(const String &,StringMap<float> &stats);
  void updateStats(GffTranscript *transcript,WigBinary &wig,
		   StringMap<float> &counts);
  void add(float x,StringMap<float> &counts,const String &word);
  void add(StringMap<float> &from,StringMap<float> &to);
  void normalize(StringMap<float> &);
  void dumpCounts(StringMap<float> &counts);
  void outputModel();
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
  : filenameRegex("([^/]+)\\.fasta"), counts(HASH_SIZE), 
    biasedCounts(HASH_SIZE)
{
  // ctor
}



void Application::main(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw BOOM::String("\n\
get-bias-model <word-length> <chunks-dir>\n\
\n\
   chunk files must be named: *.fasta *.gff *.pileup\n\
");
  wordLen=cmd.arg(0).asInt();
  String dir=cmd.arg(1);

  // Process all files in directory
  Vector<String> files;
  File::getFileList(dir,files);
  Vector<String>::iterator cur=files.begin(), end=files.end();
  for(; cur!=end ; ++cur) {
    if(filenameRegex.search(*cur)) {
      String stem=filenameRegex[1];
      processChunk(dir,stem);
    }
  }
  
  // Generate output
  normalize(counts);
  normalize(biasedCounts);
  //dumpCounts(biasedCounts);
  dumpCounts(counts);
  //outputModel();
}



void Application::processChunk(const String &dir,const String &stem)
{
  StringMap<float> countsGene(HASH_SIZE), biasedCountsGene(HASH_SIZE);

  // Get unbiased counts
  Vector<String> seqs;
  Vector<GffTranscript*> transcripts;
  getTranscriptSeqs(dir,stem,seqs,transcripts);
  int longestIndex=getLongest(seqs);
  updateStats(seqs[longestIndex],countsGene);
  normalize(countsGene);

  // Get biased counts
  String pileupFile=dir+"/"+stem+".pileup";
  WigBinary wig(pileupFile);
  updateStats(transcripts[longestIndex],wig,biasedCountsGene);
  normalize(biasedCountsGene);

  // Add normalized counts to the accumulators
  add(countsGene,counts);
  add(biasedCountsGene,biasedCounts);

  // Clean up
  Vector<GffTranscript*>::iterator cur=transcripts.begin(), 
    end=transcripts.end();
  for(; cur!=end ; ++cur) delete *cur;
}



int Application::getLongest(Vector<String> &seqs)
{
  int N=seqs.size();
  if(N==0) throw "no transcripts";
  int longest=0;
  for(int i=1 ; i<N ; ++i)
    if(seqs[i].length()>seqs[longest].length()) longest=i;
  return longest;
}



void Application::getTranscriptSeqs(const String &dir,const String &stem,
				    Vector<String> &seqs,
				    Vector<GffTranscript*> &trans)
{
  String base=dir+"/"+stem;
  String fastaFile=base+".fasta";
  String annoFile=base+".gff";
  FastaReader reader(fastaFile);
  String defline, sequence;
  reader.nextSequence(defline,sequence);
  GffReader gffReader(annoFile);
  Vector<GffTranscript*> *transcripts=gffReader.loadTranscripts();
  Vector<GffTranscript*>::iterator cur=transcripts->begin(), 
    end=transcripts->end();
  for(; cur!=end ; ++cur) {
    GffTranscript *transcript=*cur;
    transcript->loadSequence(sequence);
    String transcriptSeq=transcript->getSequence();
    //cout<<transcriptSeq<<endl<<endl; // sanity: check start/stop codons
    seqs.push_back(transcriptSeq);
    trans.push_back(transcript);
  }	
  delete transcripts;
}



void Application::updateStats(const String &S,StringMap<float> &counts)
{
  const int len=S.length();
  int end=len-wordLen;
  for(int i=0 ; i<end ; ++i) {
    String word=S.substring(i,wordLen);
    add(1,counts,word);
  }
}



void Application::add(float x,StringMap<float> &counts,const String &word)
{
  if(!counts.isDefined(word.c_str(),wordLen)) 
    counts.lookup(word.c_str(),wordLen)=x;
  else
    counts.lookup(word.c_str(),wordLen)+=x;
}



void Application::updateStats(GffTranscript *transcript,WigBinary &wig,
			      StringMap<float> &counts)
{
  Array1D<float> PSA;
  wig.buildPSA(PSA);
  int numExons=transcript->getNumExons();
  for(int i=0 ; i<numExons ; ++i) {
    BOOM::GffExon &exon=transcript->getIthExon(i);
    const String &seq=exon.getSequence();
    int end=seq.length()-wordLen;
    for(int pos=0 ; pos<end ; ++pos) {
      String word=seq.substring(pos,wordLen);
      float depth=PSA[pos+wordLen-1]-(pos>0 ? PSA[pos-1] : 0);
      depth/=wordLen;
      add(depth,counts,word);
    }
  }
}



void Application::dumpCounts(StringMap<float> &counts)
{
  StringMapIterator<float> cur=counts.begin(), end=counts.end();
  for(; cur!=end ; ++cur) {
    StringMapElem<float> elem=*cur;
    String word(elem.first,elem.len);
    float count=elem.second;
    cout<<word<<"\t"<<count<<endl;
  }
}



void Application::outputModel()
{
  StringMapIterator<float> cur=counts.begin(), end=counts.end();
  for(; cur!=end ; ++cur) {
    StringMapElem<float> elem=*cur;
    String word(elem.first,elem.len);
    float count=elem.second;
    float biasedCount=0;
    if(biasedCounts.isDefined(elem.first,elem.len))
      biasedCount=biasedCounts.lookup(elem.first,elem.len);
    float ratio=biasedCount/count;
    cout<<word<<"\t"<<ratio<<endl;
  }
}



void Application::normalize(StringMap<float> &counts)
{
  float sum=0;
  StringMapIterator<float> cur=counts.begin(), end=counts.end();
  for(; cur!=end ; ++cur) sum+=(*cur).second;
  if(sum>0)
    for(cur=counts.begin() ; cur!=end ; ++cur) (*cur).second/=sum;
}



void Application::add(StringMap<float> &from,StringMap<float> &to)
{
  StringMapIterator<float> cur=from.begin(), end=from.end();
  for(; cur!=end ; ++cur) {
    StringMapElem<float> &elem=*cur;
    if(!to.isDefined(elem.first,elem.len))
      to.lookup(elem.first,elem.len)=elem.second;
    else
      to.lookup(elem.first,elem.len)+=elem.second;
  }
}


