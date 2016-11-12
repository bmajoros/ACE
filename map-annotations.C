/****************************************************************
 map-annotation.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/CigarString.H"
#include "BOOM/GffReader.H"
#include "BOOM/FastaReader.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  int main(int argc,char *argv[]);
private:
  bool optS;
  String newSubstrate;
  CigarAlignment *alignment;
  void loadAlignment(const String &filename,bool fasta);
  void mapTranscript(const GffTranscript &,ostream &);
  void mapExon(GffExon &);
  void mapFeature(GffFeature &,ostream &);
  void projectGenes(const String &filename,ostream &);
  void projectNongenic(const String &filename,ostream &);
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
    { cerr << "STL exception caught in main:\n" << e.what() << endl; }
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}




int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"fs:");
  if(cmd.numArgs()!=3)
    throw String("\n\
map-annotations [options] <ref.gff> <cigar-file> <out.gff>\n\
  -s X : change substrate (chrom) to X\n\
  -f : cigar-file is actually a FASTA with /cigar=... on defline\n\
\n");
  const String refGffFile=cmd.arg(0);
  const String cigarFile=cmd.arg(1);
  const String outGff=cmd.arg(2);
  optS=cmd.option('s');
  if(optS) newSubstrate=cmd.optParm('s');

  // Load alignment (CIGAR)
  loadAlignment(cigarFile,cmd.option('f'));

  // Project genes, if any
  ofstream os(outGff.c_str());
  projectGenes(refGffFile,os);

  // Project non-genic elements, if any
  projectNongenic(refGffFile,os);

  delete alignment;
  return 0;
}



void Application::projectNongenic(const String &filename,ostream &os)
{
  GffReader reader(filename);
  Vector<GffFeature*> &features=*reader.loadFeatures();
  for(Vector<GffFeature*>::iterator cur=features.begin(), end=
	features.end() ; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    if(optS) feature->setSubstrate(newSubstrate);
    mapFeature(*feature,os);
    delete feature;
  }
  reader.close();
  delete &features;
}



void Application::projectGenes(const String &filename,ostream &os)
{
  GffReader reader(filename);
  Vector<GffTranscript*> &transcripts=*reader.loadTranscripts();
  for(Vector<GffTranscript*>::iterator cur=transcripts.begin(), end=
	transcripts.end() ; cur!=end ; ++cur) {
    GffTranscript *transcript=*cur;
    if(optS) transcript->setSubstrate(newSubstrate);
    mapTranscript(*transcript,os);
    delete transcript;
  }
  reader.close();
  delete &transcripts;
}



void Application::loadAlignment(const String &filename,bool isFasta)
{
  if(isFasta) {
    String defline, seq;
    FastaReader::load(filename,defline,seq);
    String id, remainder;
    FastaReader::parseDefline(defline,id,remainder);
    Map<String,String> attr;
    FastaReader::parseAttributes(remainder,attr);
    if(!attr.isDefined("cigar")) throw "No cigar element found in FASTA file";
    CigarString cigar(attr["cigar"]);
    alignment=cigar.getAlignment();
  }
  else {
    String line;
    File f(filename);
    line=f.getline();
    f.close();
    CigarString cigar(line);
    alignment=cigar.getAlignment();
  }
}



void Application::mapFeature(GffFeature &feature,ostream &os)
{
  int begin=feature.getBegin(), end=feature.getEnd();
  begin=alignment->mapApproximate(begin,DIR_NONE);
  end=alignment->mapApproximate(end,DIR_NONE);
  feature.setBegin(begin);
  feature.setEnd(end);
}



void Application::mapExon(GffExon &exon)
{
  int begin=exon.getBegin(), end=exon.getEnd();

  // These two lines map the splice sites across the
  // alignment, then use that to set exon boundaries:
  begin=alignment->mapApproximate(begin-2,DIR_NONE)+2;
  end=alignment->mapApproximate(end+1,DIR_NONE)-1;
  exon.setBegin(begin); exon.setEnd(end);
}



void Application::mapTranscript(const GffTranscript &refTrans,ostream &os)
{
  GffTranscript transcript=refTrans;
  for(Vector<GffExon*>::iterator cur=transcript.getExons(), end=
	transcript.getExonsEnd() ; cur!=end ; ++cur)
    mapExon(**cur);
  for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	transcript.getUTRend() ; cur!=end ; ++cur) {
    mapExon(**cur);
  }
  transcript.toGff(os);
}







