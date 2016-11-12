/****************************************************************
 gcf-list-samples.C
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
#include "BOOM/Pipe.H"
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/Regex.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Regex gzRegex;
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
  : gzRegex("gz$")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"rt:i:");
  if(cmd.numArgs()!=1)
    throw String("\ngcf-list-samples <in.gcf>\n\n");
  const String &gcfFilename=cmd.arg(0);

  // Process GCF file
  File *gcf=gzRegex.search(gcfFilename) ? new GunzipPipe(gcfFilename)
    : new File(gcfFilename);
  gcf->getline(); // ignore header
  while(!gcf->eof()) {
    String line=gcf->getline();
    line.trimWhitespace();
    Vector<String> fields; line.getFields(fields);
    if(fields.size()>0) cout<<fields[0]<<endl;
  }
  delete gcf;

  return 0;
}



