/**************************************************************
 ScoreAnalyzer.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "ScoreAnalyzer.H"
#include <iostream>
#include "BOOM/VectorSorter.H"
using namespace BOOM;


ScoreAnalyzer::ScoreAnalyzer(BOOM::Vector<double> &posScores,
			     BOOM::Vector<double> &negScores)
{
  // ctor

  analyze(posScores,negScores);
}



void ScoreAnalyzer::outputGraph(ostream &os)
{
  int n=precRecallGraph.size();
  for(int i=0 ; i<n ; ++i)
    {
      PrecisionRecallPair &p=precRecallGraph[i];
      os << p.precision << "\t" << p.recall << endl;
    }
}



void ScoreAnalyzer::addToList(BOOM::Vector<double> &rawScores,
			      BOOM::Vector<ScoreNode> &scores,
			      ScoreType scoreType)
{
  int n=rawScores.size();
  for(int i=0 ; i<n ; ++i)
    scores.push_back(ScoreNode(rawScores[i],scoreType));
}



void ScoreAnalyzer::analyze(BOOM::Vector<double> &posScores,
			    BOOM::Vector<double> &negScores)
{
  // Combine the pos & neg score list into a single sorted list
  BOOM::Vector<ScoreNode> allScores;
  addToList(posScores,allScores,POS_SCORE);
  addToList(negScores,allScores,NEG_SCORE);
  DirectComparator<ScoreNode> cmp;
  BOOM::VectorSorter<ScoreNode> sorter(allScores,cmp);
  sorter.sortAscendInPlace();
  int numPosScores=posScores.size();
  int numNegScores=negScores.size();
  
  // Find the overlap region
  int n=allScores.size(), leftIndex=-1, rightIndex=-1;
  int posCount=0, negCount=0;
  for(int i=0 ; i<n ; ++i)
    {
      ScoreNode &node=allScores[i];
      if(node.scoreType==POS_SCORE)
	{if(leftIndex<0) leftIndex=i;}
      else
	{
	  rightIndex=i;
	  if(leftIndex<0) ++negCount;
	}
    }

  // Compute precision-recall graph
  //for(int i=leftIndex ; i<rightIndex ; ++i)
  bestCC=0;
  bestAccuracy=0;
  negCount=0; posCount=0; //###
  for(int i=0 ; i<n ; ++i)//###
    {
      ScoreNode &thisNode=allScores[i], &nextNode=allScores[i+1];
      if(thisNode.scoreType==POS_SCORE) ++posCount;
      else ++negCount;
      double thisScore=thisNode.score, nextScore=nextNode.score;
      if(thisScore==nextScore) continue;
      double cutoff=(thisScore+nextScore)/2;
      double FN=posCount, TP=numPosScores-FN, TN=negCount, FP=numNegScores-TN;
      double precision=TP/(TP+FP), recall=TP/(TP+FN);
      double percentFP=FP/(FP+TP), percentFN=FN/(FN+TN);
      double accuracy=(TP+TN)/(TP+TN+FP+FN);
      if(accuracy>bestAccuracy) bestAccuracy=accuracy;
      double CC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
      if(CC>bestCC) bestCC=CC;
      if(TP+FP>0)
	precRecallGraph.push_back(PrecisionRecallPair(precision,recall,
						      cutoff,percentFP,
						      percentFN));
    }
  if(precRecallGraph.size()==0)
    {
      /*
	Empty precision-recall graph means nonoverlapping pos & neg
	sets, so we can achieve 100% precision and 100% recall simply
	by taking the highest neg-example score as the cutoff.
       */
      double cutoff=allScores[rightIndex].score; // rightmost negative example
      precRecallGraph.push_back(PrecisionRecallPair(1,1,cutoff,0,0));
    }
}



double ScoreAnalyzer::getBestCC() const
{
  return bestCC;
}



double ScoreAnalyzer::getBestAccuracy() const
{
  return bestAccuracy;
}



double ScoreAnalyzer::getMinRecallCutoff(double minRecall)
{
  int n=precRecallGraph.size();
  double bestPrec=precRecallGraph[0].precision;
  double bestCutoff=precRecallGraph[0].cutoff;
  for(int i=0 ; i<n ; ++i)
    {
      PrecisionRecallPair &pr=precRecallGraph[i];
      if(pr.recall>=minRecall && pr.precision>bestPrec)
	{
	  bestCutoff=pr.cutoff;
	  bestPrec=pr.precision;
	}
    }
  return bestCutoff;
}



double ScoreAnalyzer::getPercentileCutoff(double percentile)
{
  int n=precRecallGraph.size();
  int index=int(percentile*n+5/9.0);
  PrecisionRecallPair &pr=precRecallGraph[index];
  return pr.cutoff;
}



void ScoreAnalyzer::lookupCutoff(double cutoff,double &precision,
				 double &recall)
{
  int n=precRecallGraph.size();
  for(int i=0 ; i<n ; ++i)
    {
      PrecisionRecallPair p=precRecallGraph[i];
      if(cutoff<=p.cutoff)
	{
	  recall=p.recall;
	  precision=p.precision;
	  return;
	}
    }
  recall=0.0;
  precision=1.0;
}



void ScoreAnalyzer::lookupCutoff(double cutoff,double &precision,
				 double &recall,double &percentFP,
				 double &percentFN)
{
  int n=precRecallGraph.size();
  for(int i=0 ; i<n ; ++i)
    {
      PrecisionRecallPair p=precRecallGraph[i];
      if(cutoff<=p.cutoff)
	{
	  recall=p.recall;
	  precision=p.precision;
	  percentFP=p.percentFP;
	  percentFN=p.percentFN;
	  return;
	}
    }
  recall=0.0;
  precision=1.0;
}


