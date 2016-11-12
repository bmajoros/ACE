/****************************************************************
 score-content.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Sequence.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/File.H"
#include "ContentSensor.H"
using namespace std;
using namespace BOOM;

PureDnaAlphabet alphabet;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  void write(Array1D<double> &,File &,ostream &);
  bool wantText;
  String textFilename;
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
  CommandLine cmd(argc,argv,"t:");
  if(cmd.numArgs()!=3)
    throw String(
"score-content -t <*.model> <*.fasta> <out-filestem>\n\
  options:\n\
     -t <textfile> : emit text file(s) in addition to binary file\n\
");
  String modelFilename=cmd.arg(0);
  String fastaFilename=cmd.arg(1);
  String outfilestem=cmd.arg(2);
  wantText=cmd.option('t');
  if(wantText) textFilename=cmd.optParm('t');

  ContentSensor *sensor=ContentSensor::load(modelFilename);
  //sensor=sensor->reverseComplement();
  bool phased=sensor->isPhased();
  File *os[3];
  ostream *tos[3];
  if(phased) 
    for(int phase=0 ; phase<3 ; ++phase) {
      String fname=outfilestem+phase+".content";
      os[phase]=new File(fname,"w");
      String tname=textFilename+phase+".content.txt";
      if(wantText) tos[phase]=new ofstream(tname.c_str());
    }
  else {
    os[0]=new File(outfilestem+".content","w");
    String tname=textFilename+".content.txt";
    tos[0]=new ofstream(tname.c_str());
  }
  Array1D<double> arrays[3];

  FastaReader reader(fastaFilename);
  while(!reader.eof()) {
    String defline, seq;
    if(!reader.nextSequence(defline,seq)) break;
    Sequence sequence(seq,alphabet);
    int L=sequence.getLength();
    if(phased) {
      sensor->computeScores(sequence,arrays[0],arrays[1],arrays[2]);
      for(int phase=0 ; phase<3 ; ++phase) 
	write(arrays[phase],*os[phase],*tos[phase]);
    }
    else {
      sensor->computeScores(sequence,arrays[0]);
      write(arrays[0],*os[0],*tos[0]);
    }
  }
  reader.close();
  
  return 0;
}



void Application::write(Array1D<double> &array,File &os,ostream &textStream)
{
  int L=array.size();
  os<<L;
  if(wantText) textStream<<L<<endl;
  for(int i=0 ; i<L ; ++i) {
    os<<array[i];
    if(wantText) textStream<<array[i]<<endl;
  }
}


