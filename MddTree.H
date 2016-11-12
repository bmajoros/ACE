/**************************************************************
 MddTree.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_MddTree_H
#define INCL_MddTree_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Alphabet.H"
#include "SignalSensor.H"
#include "TreeNode.H"

class MddTree : public SignalSensor
{
  TreeNode *root;
  SignalSensor *pooledModel;

  void load(istream &,GarbageCollector &);
public:
  MddTree(const MddTree &,GarbageCollector &,bool revComp);
  MddTree(TreeNode *root,SignalSensor *pooledModel,GarbageCollector &,
	  Strand strand,SignalType,double cutoff,int consensusLength,
	  int consensusOffset,int contextWindowLength);
  MddTree(istream &,GarbageCollector &);
  MddTree(const BOOM::String &filename,GarbageCollector &);
  virtual ~MddTree();
  void computeSubmodels(ModelBuilder &,ModelType,SignalType,
			int consensusOffset,int consensusLength,
			int contextWindowLength);
  int getMaxDepth();
  void printOn(ostream &) const;
  void describe(ostream &);

  virtual double getLogP(const Sequence &,const BOOM::String &,int begin);
  virtual SignalSensor *reverseComplement();
  virtual bool save(const BOOM::String &filename);
  virtual bool save(ostream &os);
  virtual SignalSensor *clone() const;

  virtual void useLogOdds(SignalSensor &nullModel) {}
  virtual void useLogOdds_anonymous(ContentSensor &nullModel) {}
};

ostream &operator<<(ostream &,MddTree &);


#endif
