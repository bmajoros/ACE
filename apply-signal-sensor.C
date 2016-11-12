/****************************************************************
 apply-signal-sensor.C
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/Sequence.H"
#include "BOOM/Constants.H"
#include "BOOM/Vector.H"
#include "BOOM/VectorSorter.H"
#include "SignalSensor.H"
#include "ContentSensor.H"
#include "GarbageCollector.H"
using namespace std;
using namespace BOOM;


Alphabet alphabet=DnaAlphabet::global();


class Application
{
  Vector<double> scores;
  ContentSensor *background, *revBackground;
  int desiredCount;
  bool wantConsensusPos;
  void scan(const String &seq,const String &ID,SignalSensor *fw,
	    SignalSensor *rev,String label);
public:
  Application();
  int main(int argc,char *argv[]);
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
    CommandLine cmd(argc,argv,"c:n:b:C");
    if(cmd.numArgs()!=3) throw String(
"apply-signal-sensor [options] <sensor> <label> <*.fasta>\n\
    where options are:\n\
       -C = give consensus position, rather than context window\n\
       -c cutoff\n\
       -n #predictions\n\
       -b <infile> = use background model, score LLR's\n\
");
    String sensorFile=cmd.arg(0);
    String label=cmd.arg(1);
    String fastaFile=cmd.arg(2);
    wantConsensusPos=cmd.option('C');
    desiredCount=cmd.option('n') ? cmd.optParm('n').asInt() : -1;
    background=
      (cmd.option('b') ? ContentSensor::load(cmd.optParm('b')) : NULL);
    revBackground=(background ? background->reverseComplement() : NULL);

    // Load sensor
    GarbageIgnorer GC;
    SignalSensor *sensor=SignalSensor::load(sensorFile,GC);
    SignalSensor *revSensor=sensor->reverseComplement();
    if(cmd.option('c')) {
      if(desiredCount>=0) throw "options -c and -n are mutually exclusive";
      double cutoff=cmd.optParm('c').asDouble();
      sensor->setCutoff(cutoff);
      revSensor->setCutoff(cutoff);
    }
    if(desiredCount>=0) {
      sensor->setCutoff(NEGATIVE_INFINITY);
      revSensor->setCutoff(NEGATIVE_INFINITY);
    }

    // Process FASTA file
    FastaReader reader(fastaFile);
    String defline, seq, ID, junk;
    while(reader.nextSequence(defline,seq)) {
      reader.parseDefline(defline,ID,junk);
      scan(seq,ID,sensor,revSensor,label);
    }
    reader.close();

    // Handle the case of a known number of desired predictions
    if(desiredCount>0) { 
      DirectComparator<double> cmp;
      VectorSorter<double> sorter(scores,cmp);
      sorter.sortDescendInPlace();
      int lowerBound=NEGATIVE_INFINITY, upperBound=POSITIVE_INFINITY;
      int n=scores.size();
      double threshold;
      while(1) { // just an "if" with "breaks"
	if(n<=desiredCount) {threshold=NEGATIVE_INFINITY; break;}
	int lowerIndex=desiredCount-1, upperIndex=desiredCount;
	double a=scores[lowerIndex], b=scores[upperIndex];
	if(a!=b) {threshold=(a+b)/2; break;}
	while(lowerIndex>0 && scores[lowerIndex-1]==a) --lowerIndex;
	while(upperIndex<n-1 && scores[upperIndex+1]==a) ++upperIndex;
	int d1=desiredCount-lowerIndex, d2=upperIndex-desiredCount;
	if(d1<d2) {
	  if(lowerIndex>0) 
	    threshold=(scores[lowerIndex]+scores[lowerIndex-1])/2;
	  else threshold=0;
	  break;
	}
	if(upperIndex<n-1) {
	  threshold=(scores[upperIndex+1]+scores[upperIndex])/2;
	}
	else threshold=POSITIVE_INFINITY;
	break;
      }
      cerr<<"using threshold "<<threshold<<endl;
      desiredCount=-1;
      sensor->setCutoff(threshold);
      revSensor->setCutoff(threshold);
      FastaReader reader(fastaFile);
      String defline, seq, ID, junk;
      while(reader.nextSequence(defline,seq)) {
	reader.parseDefline(defline,ID,junk);
	scan(seq,ID,sensor,revSensor,label);
      }
    }
    return 0;
  }



void Application::scan(const String &seq,const String &ID,
		       SignalSensor *sensor,SignalSensor *revSensor,
		       String label)
{
  Sequence S(seq,alphabet);
  int sigLen=sensor->getContextWindowLength();
  int seqLen=S.getLength();
  for(int i=0 ; i<=seqLen-sigLen ; ++i) { // ### <= was <
    //SignalPtr ptr=sensor->detect(S,seq,i);
    double fgScore=sensor->getLogP(S,seq,i), score;
    double rev_fgScore=revSensor->getLogP(S,seq,i), revScore;
    if(background) {
      double bgScore=background->scoreSubsequence(S,seq,i,sigLen,0);
      double rev_bgScore=revBackground->scoreSubsequence(S,seq,i,sigLen,0);
      score=fgScore-bgScore;
      revScore=rev_fgScore-rev_bgScore;
    }
    else {
      score=fgScore;
      revScore=rev_fgScore;
    }
    if(score>=sensor->getCutoff()) {
      int begin=i+1, end=i+sigLen;
      if(wantConsensusPos) {
	begin=i+sensor->getConsensusOffset()+1;
	end=begin+sensor->getConsensusLength()-1;
      }
      if(desiredCount>=0) scores.push_back(score);
      else cout<<ID<<"\tSignalFinder\t"<<label<<"\t"<<begin<<"\t"
	       <<end<<"\t"<<score<<"\t+\t."<<endl;
    }
    if(revScore>=revSensor->getCutoff()) {
      int begin=i+1, end=i+sigLen;
      if(wantConsensusPos) {
	begin=i+revSensor->getConsensusOffset()+1;
	end=begin+revSensor->getConsensusLength()-1;
      }
      if(desiredCount>=0) scores.push_back(revScore);
      else cout<<ID<<"\tSignalFinder\t"<<label<<"\t"<<begin<<"\t"
	       <<end<<"\t"<<revScore<<"\t-\t."<<endl;
    }
  }
}



