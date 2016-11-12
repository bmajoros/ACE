/*
 refine-polya.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <string>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/FastaReader.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/Constants.H"
#include "BOOM/ProteinTrans.H"
#include "TrainingSequence.H"
#include "SignalSensor.H"
#include "GarbageCollector.H"
#include "WAM.H"
#include "WMM.H"

#define CONSENSUS_OFFSET 5
#define CONSENSUS_LENGTH 6

Alphabet alphabet; // required for linking in Partition
typedef BOOM::GffReader::TranscriptList TranscriptList;

//polya-consensus		= AATAAA|ATTAAA

class ScoredSequence : public TrainingSequence
{
  double score;
public:
  ScoredSequence(const BOOM::String &s,Alphabet &a,double score)
    : TrainingSequence(s,a), score(score) {}
  double getScore() {return score;}
};


struct ScoredSeqComp : BOOM::Comparator<TrainingSequence*>
{
  bool equal(TrainingSequence *&a,TrainingSequence *&b)
  {return 
     static_cast<ScoredSequence*>(a)->getScore() ==
     static_cast<ScoredSequence*>(b)->getScore();}
  bool greater(TrainingSequence *&a,TrainingSequence *&b)
  {return 
     static_cast<ScoredSequence*>(a)->getScore() >
     static_cast<ScoredSequence*>(b)->getScore();}
  bool less(TrainingSequence *&a,TrainingSequence *&b)
  {return 
     static_cast<ScoredSequence*>(a)->getScore() <
     static_cast<ScoredSequence*>(b)->getScore();}
};


class Application
{
  int margin;
  BOOM::Vector<BOOM::String> margins;
  SignalSensor *currentModel;
  int modelLength;
  GarbageIgnorer gc;

  void processContigs(const BOOM::String &contigsFile,
		      BOOM::Map<BOOM::String,TranscriptList*> &transcripts);
  void processForwardFeature(int featureEnd,const BOOM::String &seq);
  void processReverseFeature(int featureBegin,const BOOM::String &seq);
  void updateModel(bool echo,float keepFraction);
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
    if(cmd.numArgs()!=7)
      throw string("\n\
refine-polya <seed-model> <genes.gff> <contigs.fasta>\n\
             <outfile> <margin> <iterations> <percentile>\n\
");
    BOOM::String seedModelFile=cmd.arg(0);
    BOOM::String genesFile=cmd.arg(1);
    BOOM::String contigsFile=cmd.arg(2);
    BOOM::String outFile=cmd.arg(3);
    margin=cmd.arg(4).asInt();
    int iterations=cmd.arg(5).asInt();
    float keepFraction=cmd.arg(6).asFloat();

    alphabet=DnaAlphabet::global;

    // Load GFF
    BOOM::GffReader gffReader(genesFile);
    BOOM::Map<BOOM::String,TranscriptList*> &transcripts=
      *gffReader.loadByContig();

    // Load poly-A model
    currentModel=SignalSensor::load(seedModelFile,gc);
    modelLength=currentModel->getContextWindowLength();

    // Process each contig to extract regions following genes
    processContigs(contigsFile,transcripts);

    // Iteratively process the margins to search for best hits to 
    // current model
    for(int i=0 ; i<iterations ; ++i)
      updateModel(i==iterations-1,keepFraction);

    // Output refined model
    currentModel->save(outFile);

    delete &transcripts;
    return 0;
  }



void Application::processContigs(const BOOM::String &contigsFile,
	 BOOM::Map<BOOM::String,TranscriptList*> &transcriptMap)
{
  BOOM::FastaReader fastaReader(contigsFile);
  BOOM::String def, seq;
  while(fastaReader.nextSequence(def,seq))
    {
      BOOM::String substrate=BOOM::FastaReader::getId(def);
      if(transcriptMap.isDefined(substrate))
	{
	  TranscriptList &transcripts=*transcriptMap[substrate];
	  int n=transcripts.size();
	  for(int i=0 ; i<n ; ++i)
	    {
	      BOOM::GffTranscript *feature=transcripts[i];
	      switch(feature->getStrand())
		{
		case '+':
		  processForwardFeature(feature->getEnd(),seq);
		  break;
		case '-':
		  processReverseFeature(feature->getBegin(),seq);
		  break;
		default:
		  continue;
		}
	    }
	}
    }
}



void Application::processForwardFeature(int featureEnd,
					const BOOM::String &seq)
{
  int begin=featureEnd;
  int end=begin+margin;
  int len=seq.length();
  if(end>=len) end=len-1;
  BOOM::String subseq=seq.substr(begin,end-begin);
  margins.push_back(subseq);
}



void Application::processReverseFeature(int featureBegin,
					const BOOM::String &seq)
{
  int end=featureBegin;
  int begin=featureBegin-margin;
  if(begin<0) begin=0;
  BOOM::String subseq=
    BOOM::ProteinTrans::reverseComplement(seq.substr(begin,end-begin));
  margins.push_back(subseq);
}



void Application::updateModel(bool echo,float keepFraction)
{
  cout<<"starting with "<<margins.size()<<" margins"<<endl;
  BOOM::Vector<TrainingSequence*> windows;
  BOOM::Vector<double> scores;
  int numMargins=margins.size();
  for(int i=0 ; i<numMargins ; ++i)
    {
      BOOM::String &margin=margins[i];
      TrainingSequence marginSeq(margin,DnaAlphabet::global);
      int len=margin.length();
      int numPositions=len-modelLength;
      int bestPos;
      double bestScore=NEGATIVE_INFINITY;
      for(int pos=0 ; pos<numPositions ; ++pos)
	{
	  BOOM::String consensus=
	    margin.substring(pos+CONSENSUS_OFFSET,CONSENSUS_LENGTH);
	  if(consensus!="ATTAAA" && consensus!="AATAAA") continue;
	  double score=currentModel->getLogP(marginSeq,margin,pos);
	  if(score>=bestScore)
	    {
	      bestScore=score;
	      bestPos=pos;
	    }
	}
      if(bestScore==NEGATIVE_INFINITY) continue;



      // ### THESE TWO LINES CAUSE A LATER SEG FAULT...FIND OUT WHY!!!!
      //TrainingSequence *window=new TrainingSequence;
      //marginSeq.getSubsequence(bestPos,modelLength,*window);


      TrainingSequence *window=
	new ScoredSequence(margin.substring(bestPos,modelLength),
			   alphabet,bestScore);
      windows.push_back(window);
    }

  // Sort windows by score
  ScoredSeqComp comparator;
  BOOM::VectorSorter<TrainingSequence*> sorter(windows,comparator);
  sorter.sortDescendInPlace();

  // Keep only the top X% for retraining
  int numWindows=windows.size();
  int keepNum=int(keepFraction*numWindows);
  windows.resize(keepNum);

  // Retrain model
  cerr <<"retraining on "<<windows.size()<<" sequences..."<<endl;
  delete currentModel;
  //currentModel=new WAM(gc,windows,2,40,POLYA,0,1);
  currentModel=new WMM(gc,windows,POLYA,CONSENSUS_OFFSET,CONSENSUS_LENGTH);

  // Delete training strings
  cerr <<"deleting strings..."<<endl;
  BOOM::Vector<TrainingSequence*>::iterator cur=windows.begin(), 
    end=windows.end();
  for(; cur!=end ; ++cur) delete *cur;
}


