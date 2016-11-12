/*******************************************************************
 compute-signal-likelihoods.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*******************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/FastaReader.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/Sequence.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Constants.H"
#include "BOOM/Histogram.H"
#include "BOOM/Regex.H"
#include "SignalSensor.H"
#include "GarbageCollector.H"


Alphabet alphabet;
BOOM::Regex filestemRegex("([^/.]+).model");


class Application
{
  BOOM::Vector<double> *loadPosExamples(const BOOM::String &filename,
				      SignalSensor &);
  BOOM::Vector<double> *loadBackground(const BOOM::String &filename,
				     SignalSensor &);
  void writeHistogramFile(BOOM::Vector<double> &,const BOOM::String &filename);
  void getExtrema(BOOM::Vector<double> &scores,double &min,double &max);
  void computeRatios(const BOOM::String &filestem,Histogram<double> &posHist,
		     Histogram<double> &backgroundHist,double minScore,
		     double maxScore,int numBins);
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
    BOOM::CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=5)
      throw string(
"\ncompute-signal-likelihoods <*.model> <pos-examples.fasta> \n\
                              <neg-examples.fasta> <consensuses> <#bins>\n\
\n\
example:\n\
  compute-signal-likelihoods tag.model tag.fasta contigs.fasta TAG,TGA,TAA\n\
\n\
\n");
    BOOM::String modelFilename=cmd.arg(0);
    BOOM::String posExamples=cmd.arg(1);
    BOOM::String negExamples=cmd.arg(2);
    BOOM::String consensusArg=cmd.arg(3);
    int numBins=cmd.arg(4).asInt();
    if(!filestemRegex.search(modelFilename)) 
      throw BOOM::String("Couldn't parse filestem from filename: ")+
	modelFilename;
    BOOM::String filestem=filestemRegex[1];

    // Load the model
    alphabet=DnaAlphabet::global();
    GarbageIgnorer GC;
    SignalSensor *model=SignalSensor::load(modelFilename,GC);
    BOOM::Vector<BOOM::String> *consensuses=consensusArg.getFields(",");
    int numConsensuses=consensuses->size();
    for(int i=0 ; i<numConsensuses ; ++i)
      model->addConsensus((*consensuses)[i]);

    // Load the examples
    BOOM::Vector<double> *posScores=loadPosExamples(posExamples,*model);
    BOOM::Vector<double> *negScores=loadPosExamples(negExamples,*model);

    // Report summary stats
    BOOM::SummaryStats posStats(*posScores);
    BOOM::SummaryStats backgroundStats(*negScores);
    cout<<"positives:  "<<posStats.getMean()<<"+/-"<<posStats.getStdDev()<<
      " ("<<posStats.getMin()<<"-"<<posStats.getMax()<<")"<<endl;
    cout<<"background: "<<backgroundStats.getMean()<<"+/-"
	<<backgroundStats.getStdDev()<<" ("<<backgroundStats.getMin()
	<<"-"<<backgroundStats.getMax()<<")"<<endl;
    writeHistogramFile(*posScores,filestem+".pos-hist");
    writeHistogramFile(*negScores,filestem+".neg-hist");
    
    // Construct histograms
    double minScore=POSITIVE_INFINITY, maxScore=NEGATIVE_INFINITY;
    getExtrema(*posScores,minScore,maxScore);
    getExtrema(*negScores,minScore,maxScore);
    Histogram<double> posHist(minScore,maxScore,numBins,0.01);
    Histogram<double> backgroundHist(minScore,maxScore,numBins,1);
    posHist.addCounts(*posScores);
    backgroundHist.addCounts(*negScores);
    backgroundHist.addCounts(*posScores);

    // Compute log-likelihood ratios and write output file
    computeRatios(filestem,posHist,backgroundHist,minScore,maxScore,
      numBins);
    /*
    posHist.divideBy(backgroundHist);
    posHist.useLogs();
    BOOM::String outfile=filestem+".isp";
    posHist.save(outfile);
    */

    return 0;
  }



BOOM::Vector<double>* Application::loadPosExamples(const BOOM::String &filename,
						 SignalSensor &model)
{
  int consensusOffset=model.getConsensusOffset();
  int windowLen=model.getContextWindowLength();
  int consensusLen=model.getConsensusLength();
  BOOM::Vector<double> *scores=new BOOM::Vector<double>;
  BOOM::FastaReader reader(filename);
  BOOM::String def, seqStr;
  while(reader.nextSequence(def,seqStr))
    {
      Sequence seq(seqStr,DnaAlphabet::global());
      int seqLen=seqStr.length();
      int consensusBegin=(seqLen-consensusLen)/2;
      int windowBegin=consensusBegin-consensusOffset;
      double score=model.getLogP(seq,seqStr,windowBegin);
      if(isFinite(score)) scores->push_back(score);
    }

  return scores;
}



BOOM::Vector<double>* Application::loadBackground(const BOOM::String &filename,
						SignalSensor &model)
{
  int consensusOffset=model.getConsensusOffset();
  int windowLen=model.getContextWindowLength();
  int consensusLen=model.getConsensusLength();
  BOOM::Vector<double> *scores=new BOOM::Vector<double>;
  BOOM::FastaReader reader(filename);
  BOOM::String def, seqStr;
  while(reader.nextSequence(def,seqStr))
    {
      Sequence seq(seqStr,DnaAlphabet::global());
      int seqLen=seqStr.length();
      int lastWindowPos=seqLen-windowLen;
      for(int pos=0 ; pos<lastWindowPos ; ++pos)
	if(model.consensusOccursAt(seqStr,pos+consensusOffset))
	  {
	    double score=model.getLogP(seq,seqStr,pos);
	    if(isFinite(score)) 
	      {
		scores->push_back(score);
	      }
	  }
    }

  return scores;
}



void Application::writeHistogramFile(BOOM::Vector<double> &scores,
				     const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  BOOM::Vector<double>::iterator cur=scores.begin(), end=scores.end();
  for(; cur!=end ; ++cur)
    os<<*cur<<endl;
}



void Application::getExtrema(BOOM::Vector<double> &scores,
			     double &min,double &max)
{
  BOOM::Vector<double>::iterator cur=scores.begin(), end=scores.end();
  for(; cur!=end ; ++cur)
    {
      double score=*cur;
      if(score<min) min=score;
      if(score>max) max=score;
    }
}



void Application::computeRatios(const BOOM::String &filestem,
				Histogram<double> &posHist,
				Histogram<double> &backgroundHist,
				double minScore,
				double maxScore,int numBins)
{
  BOOM::String outfile=filestem+".isp";
  ofstream os(outfile.c_str());
  double binSize=(maxScore-minScore)/numBins;
  os<<minScore<<"\t"<<maxScore<<"\t"<<numBins<<"\t"<<binSize<<endl;
  double sum=posHist.sum();
  for(int i=0 ; i<numBins ; ++i)
    {
      double x=minScore+i*binSize;
      //double numerator=posHist.getBin(i);
      //double denominator=backgroundHist.getBin(i);
      double numerator=posHist.getBin(i)/sum;
      //(posHist.getBin(i)+backgroundHist.getBin(i));
      double denominator=1-numerator;
      double r=denominator ? numerator/denominator : 1;
      //os<<x<<"\t"<<r<<"\t"<<posHist.getBin(i)<<"\t"<<backgroundHist.getBin(i)<<endl;
      os<<x<<"\t"<<(x+binSize)<<"\t"<<log(r)<<endl;
    }
  cout<<"Ratios written into "<<outfile<<endl;
}




