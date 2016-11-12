/****************************************************************
 gff-fix-chrom-names.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/File.H"
#include "BOOM/Regex.H"
using namespace std;
using namespace BOOM;

class Application {
  Regex gffRegex;
  void fixAndReplaceFile(const String &filename);
  void fixFile(const String &filename);
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
  : gffRegex("\.gff$")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"rd");
  if(cmd.numArgs()!=1)
    throw String("\n\
gff-fix-chrom-names [options] <in.gff>\n\
\n\
  This program prepends \"chr\" to substrate name\n\
  -r = replace original file\n\
  -d = <in.gff> is a directory name, replace all files in it\n\
");
  String infile=cmd.arg(0);
  bool wantReplace=cmd.option('r');
  bool isDir=cmd.option('d');

  if(isDir) {
    Vector<String> files;
    File::getFileList(infile,files);
    for(Vector<String>::iterator cur=files.begin(), end=files.end() ;
	cur!=end ; ++cur) {
      String file=*cur;
      if(!gffRegex.search(file)) continue;
      String path=infile+"/"+file;
      fixAndReplaceFile(path);
    }
  }
  else if(wantReplace) fixAndReplaceFile(infile);
  else fixFile(infile);

  return 0;
}


void Application::fixAndReplaceFile(const String &filename)
{
  GffReader reader(filename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  ofstream os(filename.c_str());
  for(Vector<GffFeature*>::iterator cur=features->begin(), 
	end=features->end() ; cur!=end ; ++cur) {
    String gff=(*cur)->toGff();
    if(gff.length()>=3 && gff.substring(0,3)!="chr") gff=String("chr")+gff;
    os<<gff;
    delete *cur;
  }
  delete features;
}



void Application::fixFile(const String &filename)
{
  GffReader reader(filename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  for(Vector<GffFeature*>::iterator cur=features->begin(), 
	end=features->end() ; cur!=end ; ++cur) {
    String gff=(*cur)->toGff();
    if(gff.length()>=3 && gff.substring(0,3)!="chr") gff=String("chr")+gff;
    cout<<gff;
    delete *cur;
  }
  delete features;
}

