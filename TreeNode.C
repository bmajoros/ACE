/**************************************************************
 TreeNode.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "ThreePeriodicMarkovChain.H"
#include "TreeNode.H"
#include <iostream>
#include <fstream>
#include <math.h>

extern Alphabet alphabet;


TreeNode::TreeNode(const TreeNode &other)
  : left(other.left->clone()), right(other.right->clone()),
    partition(other.partition ? other.partition->clone() : NULL),
    submodel(other.submodel ? other.submodel->clone() : NULL),
    branchLogP(other.branchLogP)
{
  // copy ctor
}



TreeNode::TreeNode(Partition *partition,double logP)
  : partition(partition), left(NULL), right(NULL), submodel(NULL),
    branchLogP(logP)
{
  // ctor
}



TreeNode::TreeNode(istream &is,GarbageCollector &gc)
  : partition(NULL), left(NULL), right(NULL), submodel(NULL)
{
  // ctor

  load(is,gc);
}



TreeNode::TreeNode(const BOOM::String &filename,GarbageCollector &gc)
  : partition(NULL), left(NULL), right(NULL), submodel(NULL)
{
  // ctor

  ifstream is(filename.c_str());
  if(!is.good())
    throw BOOM::String("Error opening file ")+filename+
      "in TreeNode::TreeNode()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="MDD")
    throw BOOM::String("Attempt to load object of type ")+modelType+
      " into MDD tree";
  load(is,gc);
}



TreeNode::~TreeNode()
{
  // dtor

  delete partition;
  delete left;
  delete right;
  delete submodel;
}



BOOM::Vector<TrainingSequence*> &TreeNode::getSequences()
{
  return sequences;
}



TreeNode *TreeNode::getLeftChild()
{
  return left;
}



TreeNode *TreeNode::getRightChild()
{
  return right;
}



void TreeNode::setLeftChild(TreeNode *node)
{
  left=node;
}



void TreeNode::setRightChild(TreeNode *node)
{
  right=node;
}



Partition &TreeNode::getPartition()
{
  return *partition;
}



void TreeNode::addSequence(TrainingSequence *sequence)
{
  sequences.push_back(sequence);
}



void TreeNode::printOn(ostream &os) const
{
  printOn(os,0);
}



ostream &operator<<(ostream &os,const TreeNode &node)
{
  node.printOn(os);
  return os;
}



void TreeNode::printOn(ostream &os,int level) const
{
  for(int i=0 ; i<level ; ++i) os << "   ";
  if(partition)
    {
      os << *partition << "    P="<<exp(branchLogP)<<endl;
      
      if(left) {os<<"L:"; left->printOn(os,level+1);}
      if(right) {os << "R:";right->printOn(os,level+1);}
    }
  else 
    {
      //os << *submodel;
      os << " (" << sequences.size() << ")" << endl;
    }
}



bool TreeNode::hasPartition()
{
  return partition ? true : false;
}



void TreeNode::addSequences(BOOM::Vector<TrainingSequence*> &v)
{
  sequences.append(v);
}



bool TreeNode::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good()) throw BOOM::String("can't create file: ")+filename;
  return save(os);
}



bool TreeNode::save(ostream &os)
{
  os<<branchLogP<<endl;
  if(partition) 
    {
      os << "internal" << endl;
      partition->save(os);
      left->save(os);
      right->save(os);
    }
  else 
    {
      os << "leaf" << endl;
      submodel->save(os);
    }
  return true;
}



double TreeNode::getLogP(const Sequence &seq,const BOOM::String &str,
			 int begin)
{
  double logP=branchLogP;
  if(partition)
    switch(partition->getDirection(seq,begin))
      {
      case DIR_LEFT:
	logP+=left->getLogP(seq,str,begin);
	break;
      case DIR_RIGHT:
	logP+=right->getLogP(seq,str,begin);
	break;
      }
  else
    logP+=submodel->getLogP(seq,str,begin);
  return logP;
}



SignalSensor *TreeNode::getSubmodel()
{
  return submodel;
}




void TreeNode::getSubmodels(BOOM::Vector<SignalSensor*> &v)
{
  //throw "TreeNode::getSubmodels()";
  if(partition)
    {
      left->getSubmodels(v);
      right->getSubmodels(v);
    }
  else
    v.push_back(submodel);
}



void TreeNode::load(istream &is,GarbageCollector &gc)
{
  is>>branchLogP;
  BOOM::String nodeType;
  is >> nodeType;
  if(nodeType=="internal")
    {
      partition=new Partition(is);
      left=new TreeNode(is,gc);
      right=new TreeNode(is,gc);
    }
  else if(nodeType=="leaf")
    submodel=SignalSensor::load(is,gc);
  else throw BOOM::String("Error loading MDD tree : syntax error in model\
 file (near \"")+nodeType+"\")";
}



void TreeNode::computeSubmodels(ModelBuilder &builder,ModelType modelType,
				SignalType signalType,int consensusOffset,
				int consensusLength)
{
  if(partition)
    {
      left->computeSubmodels(builder,modelType,signalType,consensusOffset,
			     consensusLength);
      right->computeSubmodels(builder,modelType,signalType,consensusOffset,
			      consensusLength);
    }
  else
    submodel=builder.buildSignalSensor(modelType,sequences,signalType,
				       consensusOffset,consensusLength);
}



int TreeNode::getMaxDepth()
{
  int maxLeft=0, maxRight=0;
  if(left) maxLeft=left->getMaxDepth();
  if(right) maxRight=right->getMaxDepth();
  return 1+(maxLeft>maxRight ? maxLeft : maxRight);
}



TreeNode *TreeNode::reverseComplement(int sequenceLength)
{
  TreeNode *newNode;
  if(partition)
    {
      newNode=new TreeNode(partition->reverseComplement(sequenceLength),
			   branchLogP);

      // NOTE: We don't swap the left and right subtrees, because
      // we've already reverse-complemented the partition (and
      // swapping the subtrees wouldn't be right anyway)
      newNode->left=left->reverseComplement(sequenceLength);
      newNode->right=right->reverseComplement(sequenceLength);
    }
  else
    {
      newNode=new TreeNode(NULL,branchLogP);
      newNode->submodel=submodel->reverseComplement();
    }
  return newNode;
}



TreeNode *TreeNode::clone() const
{
  return new TreeNode(*this);
}


