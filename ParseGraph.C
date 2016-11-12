/****************************************************************
 ParseGraph.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "ParseGraph.H"
#ifdef EXPLICIT_GRAPHS
#include <iostream>
#include "SignalTypeProperties.H"
#include "BOOM/Vector.H"
#include "BOOM/Stack.H"
#include "BOOM/Constants.H"
#include "BOOM/Time.H"
#include "GZilla.H"
using namespace BOOM;


ParseGraph::ParseGraph(GeneZilla &genezilla)
  : genezilla(genezilla),
    vertices(new BOOM::Vector<SignalPtr >)
{
}



ParseGraph::~ParseGraph()
{
  BOOM::Vector<SignalPtr >::iterator cur=vertices->begin(), 
    end=vertices->end();
  for(; cur!=end ; ++cur) delete *cur;
  delete vertices;
}



void ParseGraph::addVertex(SignalPtr s)
{
  /*
    ### This is terribly inefficient.  It implements an insertion
        sort without even using binary search.  There is definitely a
        better way.
   */

  BOOM::Vector<SignalPtr > &vertices=*this->vertices;
  int numVertices=vertices.size();
  if(numVertices==0) vertices.push_back(s);
  else
    {
      int sPos=s->getConsensusPosition();
      for(int i=numVertices-1 ; i>=0 ; --i)
	{
	  SignalPtr t=vertices[i];
	  if(sPos>t->getConsensusPosition())
	    {
	      int sIndex=i+1;
	      vertices.push_back(NULL);
	      for(int j=numVertices ; j>sIndex ; --j)
		vertices[j]=vertices[j-1];
	      vertices[sIndex]=s;
	      return;
	    }
	}
      vertices.push_front(s);
    }
}



int ParseGraph::numVertices()
{
  return vertices->size();
}



SignalPtr ParseGraph::getIthVertex(int i)
{
  return (*vertices)[i];
}



BOOM::Stack<SignalPtr > *ParseGraph::findOptimalPath(double &parseScore)
{
  BOOM::Vector<SignalPtr> &vertices=*this->vertices;
  BOOM::Stack<SignalPtr> *path=new BOOM::Stack<SignalPtr>;

  //===================================================================
  // Use Dynamic programming to set optimal predecessors at each vertex
  //===================================================================

  BOOM::Vector<SignalPtr>::iterator cur=vertices.begin(), end=vertices.end();
  if(vertices.size()==0) return path;
  for(; cur!=end ; ++cur) {
    SignalPtr signal=*cur;
    SignalType signalType=signal->getSignalType();
    double signalScore=signal->contextWindowScore();
    Strand signalStrand=signal->getStrand();
    int defaultSignalPhase=(signalStrand==FORWARD_STRAND ? 0 : 2);
    BOOM::Set<Edge*> &edges=signal->getEdgesIn();
    int numEdges=edges.size();
    if(numEdges==0) { // left terminus... DP initialization step
      for(int i=0 ; i<3 ; ++i) {
	signal->setPredecessor(i,NULL);
	signal->getInductiveScore(i)=0;
      }
      continue;
    }
    for(int i=0 ; i<3 ; ++i) {
      signal->setPredecessor(i,NULL);
      signal->getInductiveScore(i)=NEGATIVE_INFINITY;
    }
    BOOM::Set<Edge*>::iterator eCur=edges.begin(), eEnd=edges.end();
    for(; eCur!=eEnd ; ++eCur) {
      Edge *edge=*eCur;
      SignalPtr pred=edge->getLeft();
      if(edge->isPhased()) {
	PhasedEdge *pe=static_cast<PhasedEdge*>(edge);
	if(pe->isCoding()) // exon------------------------------------
	  for(int edgePhase=0 ; edgePhase<3 ; ++edgePhase) {
	    double score=pe->getInductiveScore(edgePhase)+signalScore;
	    int signalPhase=pe->propagateForward(edgePhase);
	    if(signalType==TAG && signalPhase!=0 ||
	       signalType==NEG_ATG && signalPhase!=2)
	      continue;
	    double &bestScore=signal->getInductiveScore(signalPhase);
	    if(score>bestScore) {
	      bestScore=score;
	      signal->setPredecessor(signalPhase,pred);
	    }
	  }
	else // intron------------------------------------------------
	  for(int phase=0 ; phase<3 ; ++phase) {
	    double score=pe->getInductiveScore(phase)+signalScore;
	    double &bestScore=signal->getInductiveScore(phase);
	    if(score>bestScore) {
	      bestScore=score;
	      signal->setPredecessor(phase,pred);
	    }
	  }
      }
      else { // intergenic or UTR-----------------------------------------
	NonPhasedEdge *ne=static_cast<NonPhasedEdge*>(edge);
	double score=ne->getInductiveScore()+signalScore;
	double &bestScore=signal->getInductiveScore(defaultSignalPhase);
	if(score>bestScore) {
	  bestScore=score;
	  signal->setPredecessor(defaultSignalPhase,pred);
	}
      }
    }
  }
  
  //=======================================
  // Trace back to reconstruct optimal path
  //=======================================

  BOOM::Vector<SignalPtr >::reverse_iterator rCur=vertices.rbegin(), 
    rEnd=vertices.rend();
  double bestScore;
  SignalPtr bestTerminus=NULL;
  int bestPhase;
  for(; rCur!=rEnd ; ++rCur) { // iterate through right termini
    SignalPtr signal=*rCur;
    if(!signal->getEdgesOut().isEmpty()) break;// not a right terminus
    for(int phase=0 ; phase<3 ; ++phase) {
      double score=signal->priorInductiveScore(phase);
      if(!bestTerminus || score>bestScore) {
	bestScore=score;
	bestTerminus=signal;
	bestPhase=phase;
      }
    }
  }
  parseScore=bestScore;
  int phase=bestPhase;
  SignalPtr signal=bestTerminus;
  while(true) {
    path->push(signal);
    SignalPtr pred=signal->getPredecessor(phase);
    signal->setPredecessor((phase+1)%3,NULL);
    signal->setPredecessor((phase+2)%3,NULL);
    if(!pred) break;
    phase=genezilla.mapPhaseBack(phase,signal,pred);
    signal=pred;
  }
  return path;
}



void ParseGraph::deleteVertices(BOOM::Array1D<bool> &deleteIfTrue)
{
  BOOM::Vector<SignalPtr > &oldArray=*vertices;
  BOOM::Vector<SignalPtr > *newArray=new BOOM::Vector<SignalPtr >;
  int n=oldArray.size();
  for(int i=0 ; i<n ; ++i)
    {
      SignalPtr signal=oldArray[i];
      if(deleteIfTrue[i])
	delete signal;
      else
	newArray->push_back(signal);
    }
  
  delete vertices;
  vertices=newArray;
}



void ParseGraph::setVertexIndices()
{
  BOOM::Vector<SignalPtr > &array=*vertices;
  int n=array.size();
  for(int i=0 ; i<n ; ++i)
    array[i]->setIndex(i);
}



int ParseGraph::findSignal(SignalType type,int consensusPos) // -1 = !found
{
  // First, find any signal with the same consensusPos
  BOOM::Vector<SignalPtr> &array=*vertices;
  int arraySize=array.size();
  int min=0, max=arraySize;
  while(max-min>1) // binary search
    {
      int mid=(min+max)/2;
      SignalPtr signal=array[mid];
      int signalPos=signal->getConsensusPosition();
      if(consensusPos<signalPos) max=mid;
      else min=mid;
    }

  // Back up to the leftmost signal with this consensusPos
  while(min>0)
    {
      SignalPtr pred=array[min-1];
      if(pred->getConsensusPosition()!=consensusPos) break;
      --min;
    }
  
  // Scan through the signals having this consensusPos to find one of the
  // right type
  while(min<arraySize)
    {
      SignalPtr signal=array[min];
      if(signal->getConsensusPosition()!=consensusPos) break;
      if(signal->getSignalType()==type) return min;
      ++min;
    }

  // None found
  return -1;
}



int ParseGraph::findFirstSignalAfter(int pos) // -1 if no signals follow pos
{
  BOOM::Vector<SignalPtr> &array=*vertices;
  int arraySize=array.size();
  int min=0, max=arraySize;
  while(min<max) // binary search
    {
      int mid=(min+max)/2;
      SignalPtr signal=array[mid];
      int signalPos=signal->getConsensusPosition();
      if(signalPos<pos) min=mid+1;
      else max=mid;
    }
  if(min>=arraySize) return -1;
  //if(array[min]->getConsensusPosition()<pos) {++min;cout<<"NECESSARY!"<<endl;}
  //if(min && array[min-1]->getConsensusPosition()>=pos) throw "OK";
  return min;
}



int ParseGraph::countLeftTermini()
{
  int i=0, n=numVertices();
  for(; i<n ; ++i)
    if(!getIthVertex(i)->isLeftTerminus())
      break;
  return i;
}



int ParseGraph::countRightTermini()
{
  int n=numVertices();
  int i=n-1;
  for(; i>0 ; --i)
    if(!getIthVertex(i)->isRightTerminus())
      break;
  return n-1-i;
}



void ParseGraph::toGFF(ostream &os) const
{
  os.precision(16);
  //setPrecision(os,6);
  os<<"##gff-version 2"<<endl;
  os<<"##source-version CIA 1.0"<<endl;
  os<<"##date "<<getDateAndTime()<<endl;
  os<<"##Type SpliceGraph"<<endl;
  os<<"# stats: "<<vertices->size()<<" vertices, "
    <<countEdges()<<" edges, "
    <<genezilla.getSeqLen()<<" residues"
    <<endl;
  
  BOOM::String substrate=genezilla.getSubstrateId();

  // Emit all vertices in the graph
  BOOM::Vector<SignalPtr>::iterator cur=vertices->begin(), 
    end=vertices->end();
  for(; cur!=end ; ++cur)
    {
      // Output this vertex
      SignalPtr signal=*cur;
      int outDegree=signal->getEdgesOut().size();
      int inDegree=signal->getEdgesIn().size();
      int begin=signal->getConsensusPosition();
      int end=begin+signal->getConsensusLength();
      ++begin;
      Strand strand=signal->getStrand();
      int index=signal->getIndex();
      SignalType signalType=signal->getSignalType();
      double score=signal->contextWindowScore();
      BOOM::String source="vertex";
      String annotated=signal->isAnnotated() ? " anno=1;" : "";
      os<<substrate<<"\t"
	<<source<<"\t"
	<<signalTypeToName(signalType)<<"\t"
	<<begin<<"\t"
	<<end<<"\t"
	<<score<<"\t"
	<<strand<<"\t"
	<<"."<<"\t"
	<<"ID="<<index<<";"
	<<" in="<<inDegree<<"; out="<<outDegree<<";"
	<<annotated
	<<endl;
    }

  // Now emit the edges
  int id=0;
  for(cur=vertices->begin() ; cur!=end ; ++cur)
    {
      SignalPtr signal=*cur;
      BOOM::Set<Edge*> &edges=signal->getEdgesOut();
      BOOM::Set<Edge*>::iterator cur=edges.begin(), eend=edges.end();
      for(; cur!=eend ; ++cur)
	{
	  Edge &edge=**cur;
	  int left=edge.getLeft()->getIndex();
	  int right=edge.getRight()->getIndex();
	  int begin=edge.getFeatureBegin()+1;
	  int end=edge.getFeatureEnd();
	  ContentType type=edge.getContentType();
	  BOOM::String label=contentTypeNiceString(type);
	  Strand strand=edge.getStrand();
	  if(edge.isPhased())
	    {
	      PhasedEdge &e=(PhasedEdge&) edge;
	      for(int phase=0 ; phase<3 ; ++phase) {
		double score=e.getEdgeScore(phase);
		if(!isFinite(score)) continue;
		os<<substrate<<"\t"
		  <<"edge\t"
		  <<label<<"\t"
		  <<begin<<"\t"
		  <<end<<"\t"
		  <<score<<"\t"
		  <<strand<<"\t"
		  <<phase<<"\t"
		  <<"left="<<left<<"; "
		  <<"right="<<right<<"; "
		  <<"edgeID="<<id<<";"
		  <<endl;
	      }
	      ++id;
	    }
	  else
	    {
	      NonPhasedEdge &e=(NonPhasedEdge&) edge;
	      double score=e.getEdgeScore();
	      if(!isFinite(score)) continue;
	      os<<substrate<<"\t"
		<<"edge\t"
		<<label<<"\t"
		<<begin<<"\t"
		<<end<<"\t"
		<<score<<"\t"
		<<strand<<"\t"
		<<"."<<"\t"
		<<"left="<<left<<"; "
		<<"right="<<right<<"; "
		<<"edgeID="<<id<<";"
		<<endl;
	      ++id;
	    }
	}
    }
}



int ParseGraph::countEdges() const
{
  int n=0;
  BOOM::Vector<SignalPtr>::iterator cur=vertices->begin(), 
    end=vertices->end();
  for(; cur!=end ; ++cur)
    {
      SignalPtr signal=*cur;
      BOOM::Set<Edge*> &edges=signal->getEdgesOut();
      n+=edges.size();
    }
  return n;
}



ostream &operator<<(ostream &os,const ParseGraph &G)
{
  G.toGFF(os);
  return os;
}


/*
void ParseGraph::normalize()
{
  for(Vector<SignalPtr>::iterator cur=vertices->begin(), end=vertices->end();
      cur!=end ; ++cur) {
    SignalPtr signal=*cur;
    SignalType newType=signal->getSignalType();
    switch(newType) {
    case UTR5GT: newType=GT; break;
    case UTR5AG: newType=AG; break;
    }
    signal->setSignalType(newType);
  }
}
*/



#endif
