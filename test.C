/****************************************************************
 test.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
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


void f() { return -1; }


int Application::main(int argc,char *argv[])
{
  f();
  cout<<"ok"<<endl;
  // Process command line
  /*
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=1)
    throw String("test <infile.fasta>");
  String filename=cmd.arg(0);

  FastaReader reader(filename);
  String defline, seq, substrate, remainder;
  if(!reader.nextSequence(defline,seq)) throw filename+" : cannot read file";
  FastaReader::parseDefline(defline,substrate,remainder);
  Map<String,String> attr;
  FastaReader::parseAttributes(remainder,attr);
  cout<<attr<<endl;
  */

  /*
  Regex r("/(\\S+)=(\\S+)");
  String test("/variants=X");
  for(int i=0 ; i<100000 ; ++i) {
    test+="X ";
    int L=i+2;
    if(r.search(test)) cout<<L<<" ok"<<endl;
    else cout<<L<<" BAD"<<endl;
  }
  */

  return 0;
}

