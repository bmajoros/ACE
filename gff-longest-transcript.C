/****************************************************************
 gff-longest-transcript.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
using namespace std;
using namespace BOOM;

class Application {
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=1)
    throw String("gff-longest-transcript <in.gff>   > out.gff");
  const String infile=cmd.arg(0);

  // Load GFF file
  Vector<GffGene> &genes=*GffReader::loadGenes(infile);

  // Iterate through genes
  for(Vector<GffGene>::const_iterator cur=genes.begin(), end=genes.end() ;
      cur!=end ; ++cur) {
    const GffGene &gene=*cur;
    const int n=gene.numTranscripts();
    int longestTrans=0, longestLen=0;
    for(int i=0 ; i<n ; ++i) {
      const GffTranscript &transcript=gene.getIthTranscript(i);
      const int length=transcript.getEnd()-transcript.getBegin();
      if(length>longestLen) { longestTrans=i; longestLen=length; }
    }
    GffTranscript &transcript=gene.getIthTranscript(longestTrans);
    transcript.toGff(cout);
  }

  delete &genes;
  return 0;
}

