/****************************************************************
 TranscriptPath.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "TranscriptPath.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;

TranscriptPath::TranscriptPath()
{
  // ctor
}


void TranscriptPath::addEdge(ACEplus_Edge *e)
{
  edges.push_back(e);
  change+=e->getChange();
}



int TranscriptPath::numEdges() const
{
  return edges.size();
}



ACEplus_Edge *TranscriptPath::operator[](int i)
{
  return dynamic_cast<ACEplus_Edge*>(edges[i]);
}


double novelVertexScore(const Model &model,Vector<ACEplus_Edge*> &edges)
{
  int N=edges.size();
  if(N==0) return 0.0;
  double score=NEGATIVE_INFINITY;
  for(int i=0 ; i<N ; ++i) {
    ACEplus_Edge *edge=edges[i];
    ACEplus_Vertex *v=edge->getRight();
    if(v->isAnnotated()) continue;
    if(!isFinite(score) || v->getScore()<score) score=edge->getScore();
  }
  return score;
}



void TranscriptPath::computeScore(const Model &model)
{
  const bool INCLUDE_CONTENT=true;
  const bool INCLUDE_SIGNALS=true;
  int N=edges.size();
  if(N==0) return 0.0;
  score=INCLUDE_SIGNALS ? edges[0]->getLeft()->getScore() : 0.0;
  for(int i=0 ; i<N ; ++i) {
    ACEplus_Edge *edge=edges[i];
    if(INCLUDE_CONTENT) score+=edge->getScore();
    if(INCLUDE_SIGNALS) score+=edge->getRight()->getScore();
    if(!isFinite(edge->getScore())) cout<<*edge<<endl;
    if(!isFinite(edge->getRight()->getScore()))
      cout<<*(edge->getRight())<<endl;
  }
}



void TranscriptPath::dumpScores()
{
  int N=edges.size();
  if(N==0) return 0.0;
  score=edges[0]->getLeft()->getScore();
  cout<<"initial score="<<score<<endl;
  for(int i=0 ; i<N ; ++i) {
    ACEplus_Edge *edge=dynamic_cast<ACEplus_Edge*>(edges[i]);
    score+=edge->getScore()+edge->getRight()->getScore();
    cout<<"adding "<<edge->getScore()<<" + "<<edge->getRight()->getScore()<<endl;
  }
  cout<<"returning "<<score<<endl;
  return score;
}



ACEplus_Edge *TranscriptPath::lastEdge()
{
  return edges.empty() ? NULL : dynamic_cast<ACEplus_Edge*>(edges.getLast());
}



ACEplus_Vertex *TranscriptPath::lastVertex()
{
  ACEplus_Edge *e=lastEdge();
  return e ? dynamic_cast<ACEplus_Vertex*>(e->getRight()) : NULL;
}



TranscriptPath *TranscriptPath::clone()
{
  TranscriptPath *newPath=new TranscriptPath();
  newPath->edges=edges;
  newPath->score=score;
  newPath->change=change;
  return newPath;
}



const StructureChange &TranscriptPath::getChange() const
{
  return change;
}



GffTranscript *TranscriptPath::toTranscript(const String &transcriptID,
					    const String &geneID,
					    const String &substrate,
					    Strand strand,
					    const String &source) const
{
  GffTranscript *transcript=
    new GffTranscript(transcriptID,substrate,strand,source);
  transcript->setGeneId(geneID);
  transcript->setScore(score);
  for(Vector<ACEplus_Edge*>::const_iterator cur=edges.begin(), 
	end=edges.end() ; cur!=end ; ++cur) {
    ACEplus_Edge *edge=*cur;
    if(edge->getType()==EXON) {
      //int begin=edge->getBegin(), end=edge->getEnd();
      int begin=edge->getLeft()->getEnd();
      int end=edge->getRight()->getBegin();
      SignalType leftType=edge->getLeft()->getType();
      SignalType rightType=edge->getRight()->getType();
      if(isStartOrStop(leftType)) begin-=3;
      if(isStartOrStop(rightType)) end+=3;
      GffExon *exon=new GffExon(ET_EXON,begin,end,*transcript,false,0.0,
				false,0);
      transcript->addExon(exon);
    }
  }
  return transcript;
}



void TranscriptPath::getVertices(Vector<ACEplus_Vertex*> &into)
{
  int numEdges=edges.size();
  if(numEdges==0) return;
  into.push_back(edges[0]->getLeft());
  for(int i=0 ; i<numEdges ; ++i) {
    ACEplus_Edge *edge=edges[i];
    into.push_back(edges[i]->getRight());
  }
}



bool TranscriptPath::isFullyAnnotated() const
{
  for(Vector<ACEplus_Edge*>::const_iterator cur=edges.begin(), end=edges.end() ;
      cur!=end ; ++cur)
    if(!(*cur)->getLeft()->isAnnotated()) return false;
  return true;
}



void TranscriptPath::printOn(ostream &os) const
{
  for(Vector<ACEplus_Edge*>::const_iterator cur=edges.begin(), 
	end=edges.end() ; cur!=end ; ++cur) {
    ACEplus_Edge *edge=*cur;
    os<<*edge<<endl;
  }
}



ostream &operator<<(ostream &os,const TranscriptPath &path)
{
  path.printOn(os);
  return os;
}


