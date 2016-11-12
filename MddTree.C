/**************************************************************
 MddTree.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "MddTree.H"
#include <fstream>



MddTree::MddTree(const MddTree &other,GarbageCollector &gc,bool revcomp)
  : SignalSensor(gc,other,revcomp)
{
  // copy ctor

  if(revcomp)
    {
      int contextWindowLength=getContextWindowLength();
      int consensusLength=getConsensusLength();
      root=other.root->reverseComplement(contextWindowLength);
      pooledModel=
	(other.pooledModel ? other.pooledModel->reverseComplement() : NULL);
    }
  else
    {
      root=other.root->clone();
      pooledModel=(other.pooledModel ? other.pooledModel->clone() : NULL);
    }
}



MddTree::MddTree(TreeNode *root,SignalSensor *pooledModel,
		 GarbageCollector &gc,Strand strand,SignalType signalType,
		 double cutoff,int consensusLength,
		 int consensusOffset,int contextWindowLength)
  : root(root), 
    pooledModel(pooledModel),
    SignalSensor(gc)
{
  // ctor

  setStrand(strand);
  setSignalType(signalType);
  setSizes(consensusLength,consensusOffset,contextWindowLength);
  setCutoff(cutoff);
}



MddTree::MddTree(istream &is,GarbageCollector &gc)
  : root(NULL),
    pooledModel(NULL),
    SignalSensor(gc)
{
  // ctor

  load(is,gc);
}



MddTree::MddTree(const BOOM::String &filename,GarbageCollector &gc)
  : root(NULL), 
    pooledModel(NULL), 
    SignalSensor(gc)
{
  // ctor

  ifstream is(filename.c_str());
  if(!is.good())
    throw BOOM::String("Error loading file ")+filename+
      " in MddTree::MddTree()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="MDD")
    throw BOOM::String("Attempt to load an object of type ")+modelType
      +" into an MDD tree";
  load(is,gc);
}



MddTree::~MddTree()
{
  // dtor

  delete root;
  delete pooledModel;
}



bool MddTree::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good())
    throw BOOM::String("Error creating file ")+filename
      +" in MddTree::save()";
  return save(os);
}



bool MddTree::save(ostream &os)
{
  os << "MDD" << endl;

  os << getSignalType() << endl;
  os << getCutoff() << endl;
  os << getContextWindowLength() << '\t' 
     << getConsensusOffset() << '\t'
     << getConsensusLength() << endl;
  os << getStrand() << endl;

  root->save(os);
  if(pooledModel) 
    {
      os << 1 << endl;
      pooledModel->save(os);
    }
  else os << 0 << endl;
  return true;
}



double MddTree::getLogP(const Sequence &seq,const BOOM::String &str,
			int begin)
{
  double score;
  if(pooledModel)
    score=
      (root->getLogP(seq,str,begin)+
       pooledModel->getLogP(seq,str,begin))
      /2;
  else
    score=root->getLogP(seq,str,begin);
  return score;
}



void MddTree::load(istream &is,GarbageCollector &gc)
{
  double cutoff;
  int contextWindowLength, consensusOffset, consensusLength;
  Strand strand;
  SignalType signalType;
  BOOM::String cutoffStr;
  is >> signalType >> cutoffStr >>contextWindowLength >> consensusOffset 
     >> consensusLength >> strand;
  cutoff=cutoffStr.asDouble();
  setSignalType(signalType);
  setStrand(strand);
  setSizes(consensusLength,consensusOffset,contextWindowLength);
  setCutoff(cutoff);

  if(root) delete root;
  root=new TreeNode(is,gc);

  int usePooledModel;
  is >> usePooledModel;
  if(usePooledModel)
    pooledModel=SignalSensor::load(is,gc);
}



void MddTree::computeSubmodels(ModelBuilder &builder,ModelType modelType,
			       SignalType signalType,int consensusOffset,
			       int consensusLength,int contextWindowLength)
{
  setStrand(FORWARD_STRAND);
  setSignalType(signalType);
  setSizes(consensusLength,consensusOffset,contextWindowLength);

  root->computeSubmodels(builder,modelType,signalType,consensusOffset,
			 consensusLength);
}



void MddTree::describe(ostream &os)
{
  /*
  os << "MddTree (Maximal Dependence Decomposition):"
     << "\n\tmaximum depth: " << getMaxDepth()
     << "\n\tleaf models: ";
  TreeNode *node=root;
  while(node->hasPartition()) node=node->getLeftChild();
  node->getSubmodel()->describe(os);
  */
}



int MddTree::getMaxDepth()
{
  return root->getMaxDepth();
}



SignalSensor *MddTree::reverseComplement()
{
  return new MddTree(*this,getGC(),true);
}



SignalSensor *MddTree::clone() const
{
  return new MddTree(*this);
}



void MddTree::printOn(ostream &os) const
{
  os << *root << endl;
}



ostream &operator<<(ostream &os,MddTree &tree)
{
  tree.printOn(os);
  return os;
}

