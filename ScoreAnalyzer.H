/**************************************************************
 ScoreAnalyzer.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/

#ifndef INCL_ScoreAnalyzer_H
#define INCL_ScoreAnalyzer_H

#include "BOOM/Vector.H"
#include <iostream>

enum ScoreType
  {
    POS_SCORE, // true positives & false negatives
    NEG_SCORE  // true negatives & false positives
  };

struct ScoreNode
{
  double score;
  ScoreType scoreType;
  ScoreNode(double score,ScoreType scoreType) : 
    score(score), scoreType(scoreType) {}
  bool operator<(ScoreNode other) { return score<other.score; }
  bool operator==(ScoreNode other) { return score==other.score; }
  bool operator>(ScoreNode other) { return score>other.score; }
};

struct PrecisionRecallPair
{
  double precision, recall, cutoff, percentFP, percentFN;
  PrecisionRecallPair(double p,double r,double cutoff,double percentFP,
		      double percentFN) : 
    precision(p), recall(r), cutoff(cutoff), percentFP(percentFP),
    percentFN(percentFN) {}
};

class ScoreAnalyzer
{
  BOOM::Vector<PrecisionRecallPair> precRecallGraph;
  double bestAccuracy, bestCC;

  void analyze(BOOM::Vector<double> &posScores,BOOM::Vector<double> &negScores);
  void addToList(BOOM::Vector<double> &,BOOM::Vector<ScoreNode> &,ScoreType);
public:
  ScoreAnalyzer(BOOM::Vector<double> &posScores,BOOM::Vector<double> &negScores);
  void outputGraph(ostream &);
  double getMinRecallCutoff(double minRecall);
  double getPercentileCutoff(double percentile);//larger % -> higher cutoff
  void lookupCutoff(double cutoff,double &precision,double &recall);
  void lookupCutoff(double cutoff,double &precision,double &recall,
		    double &percentFP,double &percentFN);
  double getBestAccuracy() const;
  double getBestCC() const;
};


#endif
