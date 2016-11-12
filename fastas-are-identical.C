/****************************************************************
 fastas-are-identical.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
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
    throw String("fastas-are-identical <1.fasta> <2.fasta>");
  const String &file1=cmd.arg(0);
  const String &file2=cmd.arg(1);;

  FastaReader f1(file1), f2(file2);
  String def, seq1, seq2;
  while(f1.nextSequence(def,seq1)) {
    f2.nextSequence(def,seq2);
    if(seq1!=seq2) { cout<<"no"<<endl; return 0; }
  }
  seq2="";
  f2.nextSequence(def,seq2);
  if(seq2.length()>0) cout<<"No"<<endl;
  else cout<<"yes"<<endl;

  return 0;
}

