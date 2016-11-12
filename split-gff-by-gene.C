/****************************************************************
 split-gff-by-gene.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
using namespace std;
using namespace BOOM;



class Application {
  bool wantLongest;
  bool needToFixStop;
  void emitLongest(const GffGene &,ostream &);
  void emitAll(const GffGene &,ostream &);
  void fixStopCodons(GffGene &);
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
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"ls");
  if(cmd.numArgs()!=2)
    throw String("\nsplit-gff-by-gene [options] <in.gff> <out-dir>\n\
      -l = [L]ongest transcript only\n\
      -s = fix coordinates to include [S]top codon\n\
");
  const String inGff=cmd.arg(0);
  const String outDir=cmd.arg(1);
  wantLongest=cmd.option('l');
  needToFixStop=cmd.option('s');

  // Load GFF file
  Vector<GffGene> &genes=*GffReader::loadGenes(inGff);
  cerr<<genes.size()<<" genes loaded"<<endl;

  // Generate output
  for(Vector<GffGene>::iterator cur=genes.begin(), end=genes.end() ;
      cur!=end ; ++cur) {
    GffGene &gene=*cur;
    if(needToFixStop) fixStopCodons(gene);
    const String &id=gene.getID();
    const String outfile=outDir+'/'+id+".gff";
    ofstream os(outfile.c_str());
    if(wantLongest) emitLongest(gene,os);
    else emitAll(gene,os);
  }

  return 0;
}



void Application::emitLongest(const GffGene &gene,ostream &os)
{
  const int n=gene.numTranscripts();
  int longestL=0;
  GffTranscript *longest=NULL;
  for(int i=0 ; i<n ; ++i) {
    GffTranscript &trans=gene.getIthTranscript(i);
    const int length=trans.getEnd()-trans.getBegin();
    if(length>longestL) { longest=&trans; longestL=length; }
  }
  longest->toGff(os);
}



void Application::emitAll(const GffGene &gene,ostream &os)
{
  const int n=gene.numTranscripts();
  for(int i=0 ; i<n ; ++i) {
    GffTranscript &trans=gene.getIthTranscript(i);
    trans.toGff(os);
  }
}



void Application::fixStopCodons(GffGene &gene)
{
  const int n=gene.numTranscripts();
  for(int i=0 ; i<n ; ++i) {
    GffTranscript &trans=gene.getIthTranscript(i);
    trans.extendFinalExonBy3();
  }
}


