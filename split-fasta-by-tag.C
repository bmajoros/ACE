/****************************************************************
 split-fasta-by-tag.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/Regex.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  FastaWriter writer;
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
  if(cmd.numArgs()!=3)
    throw String("\n\split-fasta-by-tag <in.fasta> <tag> <outdir>\n\
\n\
    Deflines must have tags of the form:  /tag=value\n\
    One file per tag value will be created in directory <outdir>\n\
");
  const String infile=cmd.arg(0);
  const String tag=cmd.arg(1);
  const String outdir=cmd.arg(2);

  Regex defRegex(String("/")+tag+"="+"(\\S+)");
  FastaReader reader(infile);
  String def, seq;
  while(reader.nextSequence(def,seq)) {
    if(defRegex.search(def)) {
      String value=defRegex[1];
      String filename=outdir+"/"+value+".fasta";
      writer.appendToFasta(def,seq,filename);
    }
  }
  return 0;
}

