/**************************************************************
 TreeNode.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_TreeNode_H
#define INCL_TreeNode_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/Sequence.H"
#include "BOOM/Alphabet.H"
#include "SignalSensor.H"
#include "Partition.H"
#include "ModelBuilder.H"

class TreeNode
{
  BOOM::Vector<TrainingSequence*> sequences;
  TreeNode *left, *right;
  Partition *partition;
  SignalSensor *submodel;
  double branchLogP; // assessed when traversing from parent to this node

  void load(istream &,GarbageCollector &);
public:
  TreeNode(const TreeNode &); // does not copy over sequences
  TreeNode(Partition *,double branchLogP);
  TreeNode(istream &,GarbageCollector &);
  TreeNode(const BOOM::String &filename,GarbageCollector &);
  virtual ~TreeNode();
  bool hasPartition();
  BOOM::Vector<TrainingSequence*> &getSequences();
  TreeNode *getLeftChild();
  TreeNode *getRightChild();
  void setLeftChild(TreeNode *);
  void setRightChild(TreeNode *);
  Partition &getPartition();
  void addSequence(TrainingSequence *);
  void addSequences(BOOM::Vector<TrainingSequence*> &);
  void printOn(ostream &) const;
  void printOn(ostream &,int level) const;
  bool save(const BOOM::String &filename);
  bool save(ostream &);
  void computeSubmodels(ModelBuilder &,ModelType,SignalType,
			int consensusOffset,int consensusLength);
  double getLogP(const Sequence &,const BOOM::String &,int begin);
  SignalSensor *getSubmodel();
  void getSubmodels(BOOM::Vector<SignalSensor*> &);
  int getMaxDepth();
  TreeNode *reverseComplement(int contextWindowLength);
  TreeNode *clone() const;
};

ostream &operator<<(ostream &,const TreeNode &);

#endif
