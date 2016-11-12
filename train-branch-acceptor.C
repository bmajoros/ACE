/**************************************************************
 train-branch-acceptor.C
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
  SignalSensor *model, *negModel;
  float trainTestRatio;
  BOOM::Regex frameRegex;

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
    negModel(NULL), windowSize(5), frameRegex("/frame=(\\d+)"),
    modelType(BRANCH_ACCEPTOR)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // ====================
    // Process command line
    // ====================

    BOOM::CommandLine cmd(argc,argv,"s:gr::b:p:");
    if(cmd.numArgs()!=9)
      throw BOOM::String("\n\
train-branch-acceptor -g [-s MinSampleSize] [-r <null.model>]\n\
   <pos.fasta> <neg.fasta> <output-filestem> <%train>\n\
   <branch-point-context-length> <consensus-offset>\n\
   <min-sensitivity> <branch-order> <acceptor-order>\n\
   [-b num-boosting-iterations] [-p boost-percentile]\n\
where -g outputs a precision-recall graph (requires -e)\n\
      -r causes model to use log-odds ratios instead of log-probabilities\n\
         (the model file is optional; w/o it we use the model's rev comp)\n\
      -b <N> causes N iterations of boosting\n\
      -p <P> sets boosting percentile at P (default=0.1)\n\
      <%train> = portion of training set to use for training (vs. testing)\n\
                 (example: 0.90 = hold out 10% for testing)\n\
      <min-sensitivity> = (for example) 0.98\n\
\n");
    BOOM::String posFilename=cmd.arg(0);
    BOOM::String negFilename=cmd.arg(1);
    BOOM::String outFilestem=cmd.arg(2);
    trainTestRatio=cmd.arg(3).asFloat();
    if(cmd.option('s')) minSampleSize=cmd.optParm('s').asInt();
    alphabet=DnaAlphabet::global();
    int branchContextLength=cmd.arg(4).asInt();
    int consensusOffset=cmd.arg(5).asInt();
    float minSensitivity=cmd.arg(6).asFloat();
    int branchOrder=cmd.arg(7).asInt();
    order=cmd.arg(8).asInt();

    // =========================
    // Load training & test data
    // =========================

    cerr << "loading data..." << endl;
    BOOM::Vector<TrainingSequence*> posExamples, negExamples;
    load(posFilename,posExamples);
    int numPosExamples=posExamples.size();
    int numPosTrain=int(numPosExamples*trainTestRatio);
    int numPosTest=numPosExamples-numPosTrain;
    if(trainTestRatio<1)
      {
	if(cmd.option('b')) throw "-b requires that %train = 1.0";
	posExamples.getSubrange(0,numPosTrain-1,posTrain);
	posExamples.getSubrange(numPosTrain,numPosExamples-1,posTest);
      }
    else
      {
	cout << "Warning -- using training set for testing." << endl;
	posExamples.getSubrange(0,numPosExamples-1,posTrain);
	posExamples.getSubrange(0,numPosExamples-1,posTest);
      }
    
    cerr << "loading negative examples..." << endl;
    load(negFilename,negExamples);
    int numNegExamples=negExamples.size();
    int numNegTrain=int(numNegExamples*trainTestRatio);
    int numNegTest=numNegExamples-numNegTrain;
    if(trainTestRatio<1)
      {
	negExamples.getSubrange(0,numNegTrain-1,negTrain);
	negExamples.getSubrange(numNegTrain,numNegExamples-1,negTest);
      }
    else
      {
	negExamples.getSubrange(0,numNegExamples-1,negTrain);
	negExamples.getSubrange(0,numNegExamples-1,negTest);
      }
    
    int numIterations=1;
    if(cmd.option('b')) numIterations=cmd.optParm('b').asInt();
    float boostPercentile=
      cmd.option('p') ? cmd.optParm('p').asFloat() : 0.1;
    for(int i=0 ; i<numIterations ; ++i)
      {
	delete model;
	BOOM::Vector<double> posScores, negScores;

	// =======================
	// Build & train the model
	// =======================
	
	cerr << "training model on " << posTrain.size() 
	     << " positive examples" << endl;
	ModelBuilder modelBuilder(gc,alphabet,minSampleSize,order,
				  windowSize);
	modelBuilder.changeBranchContextLength(branchContextLength);
	modelBuilder.changeBranchOrder(branchOrder);
	model=modelBuilder.buildSignalSensor(modelType,posTrain,AG,
					     consensusOffset,
					     2);
	
	// =============================================================
	// Test the model and establish an acceptable cutoff (threshold)
	// =============================================================
	
	if(cmd.option('r')) 
	  {
	    cerr << "training negative model on " << negTrain.size()
		 << " negative examples" << endl;
	    negModel=modelBuilder.buildSignalSensor(modelType,negTrain,
						    AG,
						    consensusOffset,
						    2);
	
	    cerr << "using log odds ratios instead of probabilities" 
		 << endl;
	    BOOM::String filename=cmd.optParm('r');
	    if(filename=="")
	      model->useLogOdds(*negModel);
	    else
	      {
		ContentSensor *nullModel=ContentSensor::load(filename);
		model->useLogOdds_anonymous(*nullModel);
		delete nullModel;
	      }
	  }
	
	cout << "testing model on " << posTest.size() 
	     << " positive examples" << endl;
	evaluate(posTest,posScores);
	cout << "testing model on " << negTest.size()
	     << " negative examples" << endl;
	evaluate(negTest,negScores);
	cout<<posScores.size()<<" pos scores, "<<negScores.size()
	    <<" neg scores"<<endl;

	if(i==numIterations-1)
	  {
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
	  }

	ScoreAnalyzer analyzer(posScores,negScores);
	double cutoff=analyzer.getMinRecallCutoff(minSensitivity);
	double prec, rec, percentFP, percentFN;
	analyzer.lookupCutoff(cutoff,prec,rec,percentFP,percentFN);
	double boostCutoff=
	  analyzer.getMinRecallCutoff(1-boostPercentile);//0.9
	updateBoostCounts(posTest,posScores,boostCutoff);
	cout<<"boosting cutoff="<<boostCutoff<<endl;
	model->setCutoff(cutoff);
	cout << "************************************* Iteration #" 
	     << i << ": accuracy="
	     << analyzer.getBestAccuracy() 
	     << " precision=" << prec << " FP="
	     <<int(percentFP*100+5/9.0)<<"% FN="
	     <<int(percentFN*100+5/9.0)<<"%"
	     <<endl;

	if(i==numIterations-1)
	  {
	    // ==============================================
	    // Write out precision-recall graph, if requested
	    // ==============================================
	    
	    if(cmd.option('g'))
	      {
		cout << "writing precision-recall graph" << endl;
		ofstream os((outFilestem+".prec-recall").c_str());
		cout << "writing to file " 
		     << (outFilestem+".prec-recall").c_str()
		     <<endl;
		analyzer.outputGraph(os);
	      }
	  }
      }

    // ==============
    // Save the model
    // ==============

    cout << "saving model" << endl;
    model->save(BOOM::String(outFilestem+".model"));

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
      //double negModelScore=negModel->getLogP(*seq,str->c_str(),0);
      //###score-=negModelScore;
      delete str;
      scores.push_back(score);
    }
}



void Application::updateBoostCounts(BOOM::Vector<TrainingSequence*> &seqs,
				    BOOM::Vector<double> &scores,
				    double cutoff)
{
  int n=seqs.size(), updated=0;
  for(int i=0 ; i<n ; ++i)
    if(scores[i]<cutoff) 
      {
	seqs[i]->adjustBoostCount(1);
	++updated;
      }
  float percentUpdated=int(1000*updated/float(n)+5/9.0)/10.0;
  cout<<percentUpdated<<"% of training sequences were boosted ("
      <<updated<<"/"<<n<<")"<<endl;
}


