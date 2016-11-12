/****************************************************************
 BranchAcceptor.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "BranchAcceptor.H"
#include <iostream>
#include <fstream>
#include "BOOM/DnaAlphabet.H"


BranchAcceptor::BranchAcceptor(GarbageCollector &gc,
			       const BranchAcceptor &other,
			       bool revComp)
  : SignalSensor(gc,other,revComp),
    branchPoint(NULL),
    acceptor(NULL)
{
  // ctor

  if(!revComp)
    {
      acceptor=new WAM(gc,*other.acceptor,revComp);
      branchPoint=new WWAM(gc,*other.branchPoint,revComp);
    }
}



BranchAcceptor::~BranchAcceptor()
{
  delete acceptor;
  delete branchPoint;
}



BranchAcceptor::BranchAcceptor(GarbageCollector &gc,BOOM::String &filename)
  : SignalSensor(gc),
    branchPoint(NULL),
    acceptor(NULL)
{
  // ctor

  ifstream is(filename.c_str());
  if(!is.good()) throw BOOM::String("Error opening file ")+filename+
		   "in BranchAcceptor::BranchAcceptor()";
  BOOM::String modelType;
  is >> modelType;
  if(modelType!="BranchAcceptor") 
    throw BOOM::String("Attempt to load an object of type ")+modelType+
      "in into a BranchAcceptor";
  load(is);
}



BranchAcceptor::BranchAcceptor(GarbageCollector &gc,istream &is)
  : SignalSensor(gc),
    branchPoint(NULL),
    acceptor(NULL)
{
  // ctor

  load(is);
}



BranchAcceptor::BranchAcceptor(GarbageCollector &gc,
			       BOOM::Vector<TrainingSequence*> &seqs,
			       int branchPointOrder,
			       int acceptorOrder,
			       int branchContextLength,
			       int minSampleSize,
			       int consensusOffset, 
			       SignalType signalType)
  : SignalSensor(gc),
    branchPoint(NULL),
    acceptor(NULL)
{
  // ctor

  // Misc. initialization
  setStrand(FORWARD_STRAND);
  setSignalType(signalType);
  int sensorLength=seqs[0]->getLength();
  int acceptorContextLength=sensorLength-branchContextLength;
  setSizes(2,consensusOffset,sensorLength);   

  // Split training sequences into branch point windows and acceptor
  // windows
  BOOM::Vector<TrainingSequence*> branchPoints, acceptors;
  BOOM::Vector<TrainingSequence*>::iterator cur=seqs.begin(),
    end=seqs.end();
  for(; cur!=end ; ++cur)
    {
      TrainingSequence &S=**cur;
      TrainingSequence *branchPoint=new TrainingSequence();
      TrainingSequence *acceptor=new TrainingSequence();
      S.getSubsequence(0,branchContextLength,*branchPoint);
      S.getSubsequence(branchContextLength,acceptorContextLength,*acceptor);
      branchPoints.push_back(branchPoint);
      acceptors.push_back(acceptor);
    }

  // Train the branch point sensor & the acceptor sensor
  acceptor=new WAM(gc,acceptors,acceptorOrder,minSampleSize,
		   signalType,consensusOffset-branchContextLength,2);
  branchPoint=new WWAM(gc,branchPoints,branchPointOrder,minSampleSize,
		       signalType,0,0);

  // Delete the training subsequences
  cur=branchPoints.begin(); end=branchPoints.end();
  for(; cur!=end ; ++cur) delete *cur;
  cur=acceptors.begin(); end=acceptors.end();
  for(; cur!=end ; ++cur) delete *cur;
}



bool BranchAcceptor::save(const BOOM::String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good())
    throw BOOM::String("Error creating file ")+filename+
      "in BranchPoint::save()";
  return save(os);
}



bool BranchAcceptor::save(ostream &os)
{
  os.precision(8);
  os << "BranchAcceptor" << endl;
  os << getSignalType() << " " << getCutoff() << " " << getStrand() <<endl;
  os << getConsensusOffset() << endl;
  branchPoint->save(os);
  acceptor->save(os);
  return true;
}



void BranchAcceptor::load(istream &is)
{
  double cutoff;
  Strand strand;
  SignalType signalType;
  BOOM::String p;
  int consensusOffset;
  is >> signalType;
  is >> p; 
  cutoff=p.asDouble();
  is >> strand;
  is >> consensusOffset;
  setSignalType(signalType);
  setStrand(strand);
  setCutoff(cutoff);

  BOOM::String dummy;
  is>>dummy; // will always be "WWAM"
  branchPoint=new WWAM(getGC(),is);
  is>>dummy; // will always be "WAM"
  acceptor=new WAM(getGC(),is);

  int contextWindowLength=branchPoint->getContextWindowLength()+
    acceptor->getContextWindowLength();
  setSizes(2,consensusOffset,contextWindowLength);
}



double BranchAcceptor::getLogP(const Sequence &S,const BOOM::String &str,
			       int begin)
{
  double score;
  switch(getStrand())
    {
    case FORWARD_STRAND:
      {
	score=
	  branchPoint->getLogP(S,str,begin)+
	  acceptor->getLogP(S,str,
			    begin+branchPoint->getContextWindowLength());
      }
      break;

    case REVERSE_STRAND:
      {
	score=
	  branchPoint->getLogP(S,str,begin+
			       acceptor->getContextWindowLength())+
	  acceptor->getLogP(S,str,begin);
      }
      break;
    default: throw "bad!";
    }
  return score;
}



SignalSensor *BranchAcceptor::reverseComplement()
{
  BranchAcceptor *other=new BranchAcceptor(getGC(),*this,true);
  other->branchPoint=(WWAM*) branchPoint->reverseComplement();
  other->acceptor=(WAM*) acceptor->reverseComplement();
  return other;
}



void BranchAcceptor::useLogOdds(SignalSensor &nullModel)
{
  BranchAcceptor &other=dynamic_cast<BranchAcceptor&>(nullModel);
  branchPoint->useLogOdds(*other.branchPoint);
  acceptor->useLogOdds(*other.acceptor);
}



void BranchAcceptor::useLogOdds_anonymous(ContentSensor &nullModel)
{
  branchPoint->useLogOdds_anonymous(nullModel);
  acceptor->useLogOdds_anonymous(nullModel);
}



