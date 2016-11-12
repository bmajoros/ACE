/****************************************************************
 gff-change-substrate.C
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
  if(cmd.numArgs()!=2)
    throw String("gff-change-substrate <in&out.gff> <new-substrate>");
  const String gffFilename=cmd.arg(0);
  const String substrate=cmd.arg(1);

  // Load GFF file
  GffReader reader(gffFilename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  reader.close();

  // Modify substrates and save back to same file
  ofstream os(gffFilename.c_str());
  for(Vector<GffFeature*>::iterator cur=features->begin(), end=
	features->end() ; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    feature->setSubstrate(substrate);
    os<<feature->toGff();
  }
  
  return 0;
}

