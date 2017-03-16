/****************************************************************
 get-transitions.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/Array2D.H"
#include "SignalType.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  void updateCounts(GffTranscript &transcript,Array2D<int> &counts);
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
    throw String("get-transitions <in.gff> <outfile>");
  const String infile=cmd.arg(0);
  const String outfile=cmd.arg(1);
  
  // Initialize counts
  const int numSignals=numSignalTypes();
  Array2D<int> counts(numSignals,numSignals);
  counts.setAllTo(0);
  counts[LEFT_TERMINUS][RIGHT_TERMINUS]=1; // pseudocount

  // Read data
  GffReader reader(infile);
  TranscriptList *transcripts=reader.loadTranscripts();
  for(TranscriptList::iterator cur=transcripts->begin(), 
	end=transcripts->end() ; cur!=end ; ++cur) {
    GffTranscript *transcript=*cur;
    updateCounts(*transcript,counts);
  }

  // Compute probabilities and write to output
  ofstream os(outfile.c_str());
  Array2D<float> probs(numSignals,numSignals);
  probs.setAllTo(0.0);
  for(int from=0 ; from<numSignals ; ++from) {
    float sum=0.0;
    for(int to=0 ; to<numSignals ; ++to) sum+=counts[from][to];
    if(sum==0) continue;
    for(int to=0 ; to<numSignals ; ++to) {
      float P=float(counts[from][to])/sum;
      if(P==0.0) continue;
      os<<SignalType(from)<<" -> "<<SignalType(to)<<" : "<<P<<endl;
    }
  }


  return 0;
}



void Application::updateCounts(GffTranscript &transcript,Array2D<int> &counts)
{
  ++counts[LEFT_TERMINUS][TSS];
  ++counts[TES][RIGHT_TERMINUS];
  Vector<GffExon*> rawExons;
  transcript.getRawExons(rawExons);
  const int numExons=rawExons.size();
  if(numExons<1) INTERNAL_ERROR;
  int numInternalExons=numExons-2;
  if(numInternalExons<0) numInternalExons=0;
  const int numIntrons=numExons-1;
  const int numGT=numIntrons;
  const int numAG=numIntrons;
  if(numIntrons>0) {
    ++counts[TSS][GT];
    ++counts[AG][TES];
    counts[GT][AG]+=numIntrons;
    counts[AG][GT]+=numInternalExons;
  }
  else {
    ++counts[TSS][TES];
  }
  transcript.deleteExons(rawExons);
}




