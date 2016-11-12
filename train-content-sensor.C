/**************************************************************
 train-content-sensor.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Regex.H"
#include "BOOM/Alphabet.H"
#include "BOOM/DnaAlphabet.H"
#include "ScoreAnalyzer.H"
#include "ModelBuilder.H"
#include "ContentType.H"
#include "TrainingSequence.H"
#include "GarbageCollector.H"

Alphabet alphabet;
int frame=-1; // ### CAUTION: this is required by older code; to be removed

class Application
{
  ModelType modelType;
  int order, minSampleSize, windowSize;
  BOOM::Vector<TrainingSequence*> posTrain, negTrain, posTest, negTest;
  ContentSensor *model, *negModel;
  float trainTestRatio;
  BOOM::Regex frameRegex;

  void evaluate(BOOM::Vector<TrainingSequence*> &,BOOM::Vector<double> &scores);
  void load(BOOM::String filename,BOOM::Vector<TrainingSequence*> &);
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
    negModel(NULL), windowSize(5), frameRegex("/frame=(\\d+)")
  {
    // ctor
  }



void revCompSeqs(BOOM::Vector<TrainingSequence*> &forwardSeqs,
		      BOOM::Vector<TrainingSequence*> &revSeqs)
{
  BOOM::Vector<TrainingSequence*>::iterator cur=forwardSeqs.begin(), 
    end=forwardSeqs.end();
  for(; cur!=end ; ++cur)
    revSeqs.push_back(
      static_cast<TrainingSequence*>((*cur)->reverseComplement(alphabet)));
}



int Application::main(int argc,char *argv[])
  {
    // ====================
    // Process command line
    // ====================

    BOOM::CommandLine cmd(argc,argv,"o:s:gr::");
    if(cmd.numArgs()!=6)
      throw BOOM::String("\n\
train-content-sensor -g [-o order] [-s MinSampleSize] [-r <null.model>]\n\
   IMM|MC|3P|3PIMM <pos.fasta> <neg.fasta> <output-filestem>\n\
   <%train> <content-type>\n\
where -g outputs a precision-recall graph\n\
      -o sets the Markov order for MC, WAM, or 3P\n\
      -r causes the model to use log-odds ratios instead of probabilities\n\
         (the model file is optional; default = trained on neg.fasta)\n\
      MC = Nth-order Markov chain (1st order by default)\n\
      3P = 3-periodic variable-order (max N) Markov chain\n\
      3PIMM = 3-periodic Interpolated Markov chain (IMM)\n\
      IMM = Nth-order Interpolated Markov chain (IMM)\n\
      <%train> = portion of training set to use for training (vs. testing)\n\
                 (example: 0.90 = hold out 10% for testing)\n\
      <content-type> = \n\
         INTRON|INTERGENIC|SINGLE-EXON|INITIAL-EXON|INTERNAL-EXON|\n\
         FINAL-EXONS|THREE-PRIME-UTR|FIVE-PRIME-UTR\n\
\n");
    modelType=stringToModelType(cmd.arg(0));
    BOOM::String posFilename=cmd.arg(1);
    BOOM::String negFilename=cmd.arg(2);
    BOOM::String outFilestem=cmd.arg(3);
    trainTestRatio=cmd.arg(4).asFloat();
    if(cmd.option('o')) order=cmd.optParm('o').asInt();
    if(cmd.option('s')) minSampleSize=cmd.optParm('s').asInt();
    alphabet=DnaAlphabet::global();
    ContentType contentType=stringToContentType(cmd.arg(5));

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
	posExamples.getSubrange(0,numPosTrain-1,posTrain);
	posExamples.getSubrange(numPosTrain,numPosExamples-1,posTest);
      }
    else
      {
	posExamples.getSubrange(0,numPosExamples-1,posTrain);
	posExamples.getSubrange(0,numPosExamples-1,posTest);
      }

    // =======================
    // Build & train the model
    // =======================

    cerr << "training model on " << posTrain.size() 
	 << " positive examples from " << posFilename << endl;
    GarbageIgnorer GC;
    ModelBuilder modelBuilder(GC,alphabet,minSampleSize,order,windowSize);
    model=modelBuilder.buildContentSensor(modelType,posTrain,contentType);

    // ==============
    // Test the model
    // ==============

    BOOM::Vector<double> posScores, negScores;
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

    cerr << "training negative model on " << negTrain.size()
	 << " negative examples" << endl;
    negModel=modelBuilder.buildContentSensor(modelType,
					     negTrain,contentType);

    if(cmd.option('r')) 
      {
	cerr << "using log odds ratios instead of probabilities\n" << endl;
	BOOM::String filename=cmd.optParm('r');
	if(filename=="")
	  if(isCoding(contentType))
	    model->useLogOdds(*negModel);
	  else
	    model->useLogOdds(*model);
	else
	  {
	    ContentSensor *nullModel=ContentSensor::load(filename);
	    model->useLogOdds_anonymous(*nullModel);
	    delete nullModel;
	  }
      }

    int nPos=posTest.size(), nNeg=negTest.size();
    int numTest=min(nPos,nNeg);
    posTest.shuffle(); negTest.shuffle();
    posTest.resize(numTest); negTest.resize(numTest);

    cout << "testing model on " << posTest.size() 
	 << " positive examples" << endl;
    evaluate(posTest,posScores);
    cout << "testing model on " << negTest.size()
	 << " negative examples" << endl;
    evaluate(negTest,negScores);
    
    cout << "writing .scores file" << endl;
    ofstream os((outFilestem+".scores").c_str());
    BOOM::Vector<double>::iterator cur=posScores.begin(), 
      end=posScores.end();
    for(; cur!=end ; ++cur) os << "1\t" << *cur << endl;
    os << endl;
    cur=negScores.begin(), end=negScores.end();
    for(; cur!=end ; ++cur) os << "2\t" << *cur << endl;
    
    ScoreAnalyzer analyzer(posScores,negScores);

    // ==============================================
    // Write out precision-recall graph, if requested
    // ==============================================

    if(cmd.option('g'))
      {
	cout << "writing precision-recall graph" << endl;
	ofstream os((outFilestem+".prec-recall").c_str());
	cout << "writing to file " << (outFilestem+".prec-recall").c_str()
	     <<endl;
	analyzer.outputGraph(os);
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
      if(sequence.length()==0) continue;
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
      int len=seq->getLength();
      BOOM::String *str=seq->toString(alphabet);
      frame=seq->getPhase(); // ### this line will be going away...
      double score=
	model->scoreSubsequence(*seq,str->c_str(),0,len,seq->getPhase());

      double negModelScore=
	negModel->scoreSubsequence(*seq,str->c_str(),0,len,seq->getPhase());
      score-=negModelScore;

      delete str;
      scores.push_back(score);
    }
}


