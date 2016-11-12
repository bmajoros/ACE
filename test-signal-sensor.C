/**************************************************************
 test-signal-sensor.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Regex.H"
#include "BOOM/Alphabet.H"
#include "BOOM/DnaAlphabet.H"
#include "ScoreAnalyzer.H"
#include "ModelBuilder.H"
#include "SignalType.H"
#include "TrainingSequence.H"
#include "GarbageCollector.H"

Alphabet alphabet;
int frame=-1; // ### CAUTION: this is required by older code; to be removed

class Application
{
  GarbageCollector gc;
  ModelType modelType;
  int order, minSampleSize, windowSize;
  BOOM::Vector<TrainingSequence*> posTrain, negTrain, posTest, negTest;
  SignalSensor *model;
  float trainTestRatio;
  BOOM::Regex frameRegex;
  BOOM::Vector<double> posScores, negScores;

  void evaluate(BOOM::Vector<TrainingSequence*> &,BOOM::Vector<double> &scores);
  void load(BOOM::String filename,BOOM::Vector<TrainingSequence*> &);
  void updateBoostCounts(BOOM::Vector<TrainingSequence*> &seqs,
			 BOOM::Vector<double> &scores,double cutoff);
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
    catch(const BOOM::String &msg)
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
  : order(5), minSampleSize(175), model(NULL),
    windowSize(5), frameRegex("/frame=(\\d+)")
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // ====================
    // Process command line
    // ====================

    BOOM::CommandLine cmd(argc,argv,"o:s:gr::b:");
    if(cmd.numArgs()!=4)
      throw BOOM::String("\ntest-signal-sensor <*.model> <pos.fasta> <neg.fasta> <outfile>\n");
    BOOM::String modelFilename=cmd.arg(0);
    BOOM::String posFilename=cmd.arg(1);
    BOOM::String negFilename=cmd.arg(2);
    BOOM::String outFilestem=cmd.arg(3);
    alphabet=DnaAlphabet::global;

    // =========================
    // Load model
    // =========================
    model=SignalSensor::load(modelFilename,gc);

    // =========================
    // Load test data
    // =========================
    cerr << "loading data..." << endl;
    BOOM::Vector<TrainingSequence*> posExamples, negExamples;
    load(posFilename,posExamples);
    int numPosExamples=posExamples.size();
    load(negFilename,negExamples);
    int numNegExamples=negExamples.size();

    // =========================
    // Apply model
    // =========================
    cout << "testing model on " << posExamples.size() 
	 << " positive examples" << endl;
    evaluate(posExamples,posScores);
    cout << "testing model on " << negExamples.size()
	 << " negative examples" << endl;
    evaluate(negExamples,negScores);
    cout<<posScores.size()<<" pos scores, "<<negScores.size()<<" neg scores"<<endl;

    // =========================
    // Write out .scores file
    // =========================
    cout << "writing .scores file" << endl;
    ofstream os((outFilestem+".scores").c_str());
    BOOM::Vector<double>::iterator cur=posScores.begin(), 
      end=posScores.end();
    
    double minScore=*cur;
    for(; cur!=end ; ++cur) 
      if(!isinf(*cur) && *cur<minScore) minScore=*cur;
    cur=negScores.begin(); end=negScores.end();
    for(; cur!=end ; ++cur) 
      if(!isinf(*cur) && *cur<minScore) minScore=*cur;
    
    cur=posScores.begin(); end=posScores.end();	
    for(; cur!=end ; ++cur) 
      os << "1\t" << (isinf(*cur)?minScore:*cur) << endl;
    os << endl;
    cur=negScores.begin(); end=negScores.end();
    for(; cur!=end ; ++cur) 
      os << "2\t" << (isinf(*cur)?minScore:*cur) << endl;

    // =========================
    // Report accuracy
    // =========================
    ScoreAnalyzer analyzer(posScores,negScores);
    double cutoff=model->getCutoff();//analyzer.getMinRecallCutoff(minSensitivity);
    double prec, rec, percentFP, percentFN;
    analyzer.lookupCutoff(cutoff,prec,rec,percentFP,percentFN);
    updateBoostCounts(posTest,posScores,analyzer.getPercentileCutoff(0.8));
    cout<<"cutoff="<<cutoff<<endl;
    model->setCutoff(cutoff);
    cout << "accuracy="
	 << analyzer.getBestAccuracy() 
	 << " precision=" << prec << " FP="
	 <<int(percentFP*100+0.5)<<"% FN="
	 <<int(percentFN*100+0.5)<<"%"
	 <<endl;
    
    // ==============================================
    // Write out precision-recall graph
    // ==============================================
    cout << "writing precision-recall graph" << endl;
    ofstream os2((outFilestem+".prec-recall").c_str());
    cout << "writing to file " 
	 << (outFilestem+".prec-recall").c_str()
	 <<endl;
    analyzer.outputGraph(os2);
    
    cout << "done." << endl;
    return 0;
  }



void Application::load(BOOM::String filename,BOOM::Vector<TrainingSequence*> &v)
{
  BOOM::FastaReader reader(filename);
  BOOM::String defline, sequence;
  while(reader.nextSequence(defline,sequence))
    {
      if(frameRegex.search(defline))
	frame=frameRegex[1].asInt();
      else
	frame=0;
      TrainingSequence *s=new TrainingSequence(sequence,alphabet);
      s->setPhase(frame);
      v.push_back(s);
    }
}



void Application::evaluate(BOOM::Vector<TrainingSequence*> &sequences,
			   BOOM::Vector<double> &scores)
{
  int n=sequences.size();
  for(int i=0 ; i<n ; ++i)
    {
      TrainingSequence *seq=sequences[i];
      BOOM::String *str=seq->toString(alphabet);
      frame=sequences[i]->getPhase();
      double score=model->getLogP(*seq,str->c_str(),0);
      delete str;
      scores.push_back(score);
    }
}



void Application::updateBoostCounts(BOOM::Vector<TrainingSequence*> &seqs,
				    BOOM::Vector<double> &scores,
				    double cutoff)
{
  int n=seqs.size();
  for(int i=0 ; i<n ; ++i)
    if(scores[i]<cutoff) 
      seqs[i]->adjustBoostCount(1);
}


