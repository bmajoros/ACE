/**************************************************************
 mdd.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/CommandLine.H"
#include "BOOM/Stacktrace.H"
#include "BOOM/Chi2IndepTest.H"
#include "BOOM/FastaReader.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Constants.H"
#include "BOOM/Sequence.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/Environment.H"
#include "ModelBuilder.H"
#include "MddTree.H"
#include "Partition.H"
#include "TreeNode.H"
#include "ScoreAnalyzer.H"
#include "GarbageCollector.H"
#include "WMM.H"

int frame=-1;
Alphabet alphabet;

enum InferenceType
  {
    RELATIVE_ENTROPY,
    MAX_CHI
  };

inline float ln(float p) {return log(p)/log(2.0);}

class Application
{
  InferenceType inferenceType;
  int numTrain;
  SignalType signalType;
  SignalSensor *windowModel;
  int WL; // size of context window
  int consensusLength; // length of ATG/TAG/TGA/TAA/GT/AG (ie, usu. 2 or 3)
  int consensusOffset;
  ModelType windowModelType;
  int maxTreeDepth;
  int minSampleSize, order; // for Markov chains
  int windowSize; // for WWAM's only! (caution: a different type of window!)
  BOOM::Vector<TrainingSequence*> examples, negativeExamples;
  BOOM::Vector<TrainingSequence*> *trainPos, *trainNeg, *testPos, *testNeg;
  int seqLen;
  BOOM::Vector<Partition*> allPartitions;
  BOOM::Chi2Table chiSquaredTable;
  float trainTestRatio;
  BOOM::String outFilestem, negFastaFilename;
  ModelBuilder *modelBuilder;
  double minRecall, cutoff;

  void loadTrainingData(BOOM::String fastaFilename);
  void processParms(BOOM::String);
  void evaluate(SignalSensor &);
  void evaluate(SignalSensor &,
		BOOM::Vector<TrainingSequence *> &,
		BOOM::Vector<double> &scores);
  float assessPosition(Partition &,BOOM::Vector<TrainingSequence*> &,
		       int consensusPos,int otherPos);
  float assessPartition(Partition &,BOOM::Vector<TrainingSequence*> &);
  void generateAllPartitions();
  TreeNode *buildTree(BOOM::Vector<TrainingSequence*> &,int level=0,
		      double branchLogP=0);
  Partition *pickPartitionRelEnt(BOOM::Vector<TrainingSequence*> &);
  Partition *pickPartition(BOOM::Vector<TrainingSequence*> &);
  void trainNullModel(BOOM::String filename);
  void pickCutoffFromTrainSet(SignalSensor &);
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
	//throw BOOM::Stacktrace(msg.c_str());
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
  : WL(0), 
    modelBuilder(NULL), 
    cutoff(NEGATIVE_INFINITY),
    consensusLength(-1), 
    order(5), 
    minSampleSize(175)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    BOOM::CommandLine cmd(argc,argv,"d:s:o:w:p");
    if(cmd.numArgs()!=9)
      throw BOOM::String(
"mdd <consensus-offset> <consensus-length> \n"
"    <model-type> <signal-type> <pos.fasta> <neg.fasta> <out-filestem> \n"
"    <min-recall> <%train> [-d max-depth][-s #][-o #][-w #] \n"
"\n"
"     -d = maximum MDD tree depth\n"
"     -s = min sample size for Markov chains\n"
"     -o = max order for Markov chains\n"
"     -w = \"windowing\" size for WWAM\n"
"     -p = incorporate a pooled model\n"
"     model types: WMM,WAM,WWAM\n"
"     <%train> = portion of training set used for training (vs. testing)\n"
);
    consensusOffset=cmd.arg(0).asInt();
    consensusLength=cmd.arg(1).asInt();
    windowModelType=stringToModelType(cmd.arg(2));
    signalType=stringToSignalType(cmd.arg(3));
    BOOM::String fastaFilename=cmd.arg(4);
    negFastaFilename=cmd.arg(5);
    outFilestem=cmd.arg(6);
    minRecall=cmd.arg(7).asFloat();
    trainTestRatio=cmd.arg(8).asFloat();
    if(cmd.option('d')) maxTreeDepth=cmd.optParm('d').asInt();
    else maxTreeDepth=-1;
    if(cmd.option('o')) order=cmd.optParm('o').asInt();
    else cerr << "WARNING: using default order=" << order << endl;
    if(cmd.option('w')) windowSize=cmd.optParm('w').asInt();
    if(cmd.option('s')) minSampleSize=cmd.optParm('s').asInt();
    else cerr << "WARNING: using default minSampleSize=" << minSampleSize
	      << endl;
    inferenceType=MAX_CHI;
    bool usePooledModel=cmd.option('p');

    alphabet=DnaAlphabet::global();

    generateAllPartitions();
    GarbageCollector garbageCollector;
    modelBuilder=new ModelBuilder(garbageCollector,alphabet,minSampleSize,
				  order,windowSize);

    // Load training data & train context models
    cerr << "loading training data..." << endl;
    loadTrainingData(fastaFilename);

    // Train a pooled model
    SignalSensor *pooledModel=usePooledModel ?
      modelBuilder->buildSignalSensor(windowModelType,*trainPos,
				      signalType,consensusOffset,
				      consensusLength)
      : NULL;

    // Build tree
    cerr << "building tree from " << trainPos->size() 
         << " example sequences" << endl;
    TreeNode *root=buildTree(*trainPos);
    if(!root)
    {
      cerr << "WARNING: Could not apply MDD due to inadequate sample size"
	   << endl;
      root=new TreeNode(NULL,0);
      root->addSequences(*trainPos);
    }
    MddTree *tree=new MddTree(root,pooledModel,garbageCollector,
			      FORWARD_STRAND,signalType,0,consensusLength,
			      consensusOffset,WL);
    windowModel=tree;
    
    // Build leaf models
    cerr << "building leaf models..." << endl;
    tree->computeSubmodels(*modelBuilder,windowModelType,signalType,
			   consensusOffset,consensusLength,WL);

    // Train a null model, for comparison
    cerr << "training null model" << endl;
    trainNullModel(negFastaFilename);

    // Evaluate accuracy
    if(trainTestRatio<1.0)
    {
      cerr << "evaluating accuracy..." << endl;
      evaluate(*tree);
    }
    else pickCutoffFromTrainSet(*tree);
    tree->setCutoff(cutoff);

    // Save models
    cerr << "saving model" << endl;
    tree->save(outFilestem+".model");

    int depth=tree->getMaxDepth();
    cout << "tree depth=" << depth << endl;
    tree->printOn(cout);
    cout<<endl;

    delete modelBuilder;
    return 0;
  }



void Application::loadTrainingData(BOOM::String fastaFilename)
{
  BOOM::FastaReader reader(fastaFilename);
  BOOM::String defline, sequence;

  while(reader.nextSequence(defline,sequence))
    examples.push_back(new TrainingSequence(sequence,alphabet));

  seqLen=examples[0]->getLength();
  WL=seqLen;
  int numExamples=examples.size();
  numTrain=int(numExamples*trainTestRatio);

  trainPos=examples.getSubrange(0,numTrain-1);
  testPos=examples.getSubrange(numTrain,numExamples-1);
}



TreeNode *Application::buildTree(BOOM::Vector<TrainingSequence*> &sequences,
				 int level,double branchLogP)
{
  if(sequences.size()<minSampleSize || 
     level>=maxTreeDepth && maxTreeDepth>=0) 
    return NULL;

  Partition *partition=NULL;
  switch(inferenceType)
    {
    case RELATIVE_ENTROPY:
      partition=pickPartitionRelEnt(sequences);
      break;
    case MAX_CHI:
      partition=pickPartition(sequences);
      break;
    default:
      throw "ERROR!";
    }
  if(!partition) return NULL;
  TreeNode *node=new TreeNode(partition,branchLogP);
  node->addSequences(sequences);

  BOOM::Vector<TrainingSequence*> leftSeq;
  BOOM::Vector<TrainingSequence*> rightSeq;
  double leftP=partition->split(sequences,leftSeq,rightSeq);
  double rightP=1-leftP;
  if(leftSeq.size()<minSampleSize || rightSeq.size()<minSampleSize)
    {
      delete node;
      return NULL;
    }

  TreeNode *leftChild=buildTree(leftSeq,level+1,log(leftP));
  if(!leftChild) 
    {
      leftChild=new TreeNode(NULL,log(leftP));
      leftChild->addSequences(leftSeq);
    }
  TreeNode *rightChild=buildTree(rightSeq,level+1,log(rightP));
  if(!rightChild) 
    {
      rightChild=new TreeNode(NULL,log(rightP));
      rightChild->addSequences(rightSeq);
    }
  node->setLeftChild(leftChild);
  node->setRightChild(rightChild);

  return node;
}



void Application::generateAllPartitions()
{
  int numBits=alphabet.getNumElements();
  int maxCounter=(int)pow(2.0,numBits)-1;
  for(unsigned int counter=1 ; counter<maxCounter ; ++counter)
    {
      Partition *partition=new Partition;
      unsigned int mask=1;
      for(unsigned int bit=0 ; bit<numBits ; ++bit)
	{
	  if(counter & mask)
	    partition->addSymbol(Symbol(bit));
	  mask <<= 1;
	}
      allPartitions.push_back(partition);
    }
}



Partition *Application::pickPartitionRelEnt(BOOM::Vector<TrainingSequence*> &
					    sequences)
{
  float bestScore=0;
  Partition *bestPartition=NULL;
  int bestPos;

  float numSeq=sequences.size();
  int portion=int(numSeq/3);
  int numPartitions=allPartitions.size();
  for(int pos=0 ; pos<WL ; ++pos)
    {
      for(int i=0 ; i<numPartitions ; ++i)
	{
	  Partition *partition=allPartitions[i];
	  partition->setIndex(pos);
	  BOOM::Vector<TrainingSequence*> leftSeq, rightSeq;
	  partition->split(sequences,leftSeq,rightSeq);
	  float leftPortion=leftSeq.size()/numSeq;
	  float rightPortion=rightSeq.size()/numSeq;
	  if(leftPortion==0.0 || rightPortion==0.0) continue;
	  float splitEntropy=
	    -1*(leftPortion*ln(leftPortion)+rightPortion*ln(rightPortion));
	  //if(leftSeq.size()<portion || rightSeq.size()<portion) continue;
	  GarbageIgnorer gc;
	  SignalType st=NO_SIGNAL_TYPE;
	  WMM leftMatrix(gc,leftSeq,st,consensusOffset,consensusLength);
	  WMM rightMatrix(gc,rightSeq,st,consensusOffset,consensusLength);
	  float score=splitEntropy*leftMatrix.divergence(rightMatrix);
	  if(!bestPartition || score>bestScore)
	    {
	      bestPartition=partition;
	      bestScore=score;
	      bestPos=pos;
	    }
	}
    }
  if(bestPartition==NULL || bestScore==0) return NULL;
  return new Partition(*bestPartition,bestPos);  
}



Partition *Application::pickPartition(BOOM::Vector<TrainingSequence*> 
				      &sequences)
{
  float bestScore=0;
  Partition *bestPartition=NULL;
  int bestPos;

  //int portion=sequences.size()/3; //### NOT IN GENSCAN!
  int numPartitions=allPartitions.size();
  for(int pos=0 ; pos<WL ; ++pos)
    {
      for(int i=0 ; i<numPartitions ; ++i)
	{
	  Partition *partition=allPartitions[i];
	  partition->setIndex(pos);

	  BOOM::Vector<TrainingSequence*> leftSeq, rightSeq;
	  partition->split(sequences,leftSeq,rightSeq);
	  //if(leftSeq.size()<portion || rightSeq.size()<portion) 
	  //  continue; // ### NOT IN GENSCAN!

	  float score=assessPartition(*partition,sequences);
	  if(!bestPartition || score>bestScore)
	    {
	      bestPartition=partition;
	      bestScore=score;
	      bestPos=pos;
	    }
	}
    }
  if(bestPartition==NULL || bestScore==0) return NULL;
  return new Partition(*bestPartition,bestPos);
}



float Application::assessPartition(Partition &partition,
				  BOOM::Vector<TrainingSequence*> &sequences)
{
  int df=alphabet.getNumElements()-2;
  float totalChiSquared=0;
  int numSignificant=0;
  int index=partition.getIndex();
  for(int pos=0 ; pos<WL ; ++pos)
    {
      if(pos==index) continue;
      float chiSquared=assessPosition(partition,sequences,index,pos);
      totalChiSquared+=chiSquared;
      float P=chiSquaredTable.lookupP(df,chiSquared);
      if(P<=0.05) ++numSignificant;
    }
  if(numSignificant==0) return 0;
  return totalChiSquared;
}



float Application::assessPosition(Partition &partition,
				  BOOM::Vector<TrainingSequence*> &sequences,
				  int consensusPos,int otherPos)
{
  BOOM::IntArray2D contingencyTable(alphabet.getNumElements()-1,2);
  contingencyTable.setAllTo(0);
  Symbol N=alphabet.lookup('N');

  int n=sequences.size();
  for(int i=0 ; i<n ; ++i)
    {
      Sequence &seq=*sequences[i];
      Symbol consensusSymbol=seq[consensusPos];
      Symbol otherSymbol=seq[otherPos];
      if(otherSymbol==N) continue;
      if(otherSymbol>N) otherSymbol=Symbol(otherSymbol-1);
      int consensusIndicator=
	partition.isLeftSymbol(consensusSymbol) ? 1 : 0;
      int symbolCode=otherSymbol;
      ++contingencyTable[symbolCode][consensusIndicator];
    }
  BOOM::Chi2IndepTest test(contingencyTable,chiSquaredTable);
  return test.getChiSquared();
}



void Application::trainNullModel(BOOM::String filename)
{
  // Load negative examples
  BOOM::Vector<TrainingSequence*> negWindows, negLefts, negRights;
  BOOM::FastaReader reader(filename);
  BOOM::String defline, sequence;
  while(reader.nextSequence(defline,sequence))
    negativeExamples.push_back(new TrainingSequence(sequence,alphabet));
  int numNeg=negativeExamples.size();
  int numTrain=int(numNeg*trainTestRatio);

  // Split into train/test sets
  trainNeg=negativeExamples.getSubrange(0,numTrain-1);
  testNeg=negativeExamples.getSubrange(numTrain,negativeExamples.size()-1);

  // Don't train neg models on more data than positive models (saves time)
  int trainPosSize=trainPos->size();
  int testPosSize=testPos->size();
  if(trainNeg->size()>trainPosSize)
    trainNeg->resize(trainPosSize);
  if(testNeg->size()>testPosSize)
    testNeg->resize(testPosSize);
}



void Application::pickCutoffFromTrainSet(SignalSensor &signalModel)
{
  BOOM::Vector<double> posScores, negScores;
  evaluate(signalModel,negativeExamples,negScores);
  evaluate(signalModel,examples,posScores);

  // Analyze results (precision & recall)
  ScoreAnalyzer scoreAnalyzer(posScores,negScores);
  cutoff=scoreAnalyzer.getMinRecallCutoff(minRecall);
  cout << "cutoff=" << cutoff << endl;
}


void Application::evaluate(SignalSensor &signalModel)
{
  // Perform evaluation
  BOOM::Vector<double> posScores, negScores;
  cerr << "evaluating on " << testNeg->size() << " negative test cases"
       << endl;
  evaluate(signalModel,/*negativeExamples*/*testNeg,negScores);
  cerr << "evaluating on " << testPos->size() << " positive test cases"
       << endl;
  evaluate(signalModel,/*examples*/*testPos,posScores);

  // Analyze results (precision & recall)
  ofstream ospr((outFilestem+".prec-recall").c_str());
  ScoreAnalyzer scoreAnalyzer(posScores,negScores);
  scoreAnalyzer.outputGraph(ospr);
  ospr.close();
  cutoff=scoreAnalyzer.getMinRecallCutoff(minRecall);
  double precision, recall, percentFP, percentFN;
  scoreAnalyzer.lookupCutoff(cutoff,precision,recall,percentFP,percentFN);
  cout << "cutoff=" << cutoff << ", precision="
       << precision << ", recall=" << recall << ", FP="
       << int(percentFP*100+5/9.0) << "%, FN="
       << int(percentFN*100+5/9.0) << "%" << endl;
  signalModel.setCutoff(cutoff);

  // Output raw scores
  ofstream os((outFilestem+".scores").c_str());
  BOOM::Vector<double>::iterator cur=posScores.begin(), 
    end=posScores.end();
  for(; cur!=end ; ++cur) os << "1\t" << *cur << endl;
  os << endl;
  cur=negScores.begin(), end=negScores.end();
  for(; cur!=end ; ++cur) os << "2\t" << *cur << endl;
  os.close();
}



void Application::evaluate(SignalSensor &signalModel,
			   BOOM::Vector<TrainingSequence *> &examples,
			   BOOM::Vector<double> &scores)
{
  int numExamples=examples.size();
  for(int j=0 ; j<numExamples ; ++j)
    {
      TrainingSequence &sequence=*examples[j];
      BOOM::String *str=sequence.toString(alphabet);
      double score=
	signalModel.getLogP(sequence,str->c_str(),0);
      delete str;
      scores.push_back(score);
    }
}

