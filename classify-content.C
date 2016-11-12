/****************************************************************
 classify-content.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/FastaReader.H"
#include "ContentSensor.H"
using namespace std;
using namespace BOOM;

PureDnaAlphabet alphabet;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  bool useLLR;
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"n");
  if(cmd.numArgs()!=3)
    throw String("\n\
classify-content <pos.model> <neg.model> <test.fasta>\n\
   -n = ignore negative model\n\
");
  String posFilename=cmd.arg(0);
  String negFilename=cmd.arg(1);
  String testFilename=cmd.arg(2);
  useLLR=!cmd.option('n');

  ContentSensor *posModel=ContentSensor::load(posFilename);
  ContentSensor *negModel=useLLR ? ContentSensor::load(negFilename) : NULL;
  FastaReader reader(testFilename,alphabet);
  String defline, seqstr;
  while(reader.nextSequence(defline,seqstr)) {
    Sequence seq(seqstr,alphabet);
    int L=seq.getLength();
    double posScore=
      posModel->scoreSubsequence(seq,seqstr,0,L,0);
    if(useLLR) {
      double negScore=
	negModel->scoreSubsequence(seq,seqstr,0,L,0);
      double ratio=exp(posScore-negScore);
      cout<<ratio<<endl;
    }
    else cout<<posScore/L<<endl;
  }
  
  return 0;
}

