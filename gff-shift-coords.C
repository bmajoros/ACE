/****************************************************************
 gff-shift-coords.C
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
  CommandLine cmd(argc,argv,"r");
  if(cmd.numArgs()!=2)
    throw String("\n\
gff-shift-coords [OPTIONS] <in.gff> -- <delta>\n\
   The -- allows you to use a negative value for <delta>\n\
   -r = replace original file (otherwise, dump to STDOUT)\n\
");
  const String &inGff=cmd.arg(0);
  const int delta=cmd.arg(1).asInt();
  bool wantReplace=cmd.option('r');

  // Load GFF
  GffReader reader(inGff);
  Vector<GffFeature*> *features=reader.loadFeatures();

  // Fix coords
  for(Vector<GffFeature*>::iterator cur=features->begin(), 
	end=features->end() ; cur!=end ; ++cur)
    (*cur)->shiftCoords(delta);

  // Output results
  ofstream os;
  if(wantReplace) os.open(inGff.c_str());
  for(Vector<GffFeature*>::iterator cur=features->begin(), 
	end=features->end() ; cur!=end ; ++cur) {
    String gff=(*cur)->toGff();
    (wantReplace ? os : cout) << gff;
    delete *cur;
  }
  delete features;

  return 0;
}

