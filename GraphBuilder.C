/****************************************************************
 GraphBuilder.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "GraphBuilder.H"
#include "BOOM/Constants.H"
#include "BOOM/Interval.H"
#include "ACEplus_Edge.H"
#include "ACEplus_Vertex.H"
#include "SignalPrinter.H"
#include "BOOM/Stack.H"
using namespace std;
using namespace BOOM;

/* NOTES:
   1. For vertices, the begin and end positions delimit the consensus,
      not the full scoring window ("context window").
   2. For edges, the begin and end coordinates delimit the boundaries of the
      context windows of the incident vertices, not their consensuses.
   3. Thus, the coordinates of a vertex are not the same as the coordinates
      of its incident edges.
 */

#define SANITY_CHECKS


GraphBuilder::GraphBuilder(const GffTranscript &projected,
			   const TranscriptSignals &signals,
			   Model &model,
			   const Sequence &refSeq,const String &refSeqStr,
			   const Sequence &altSeq,const String &altSeqStr,
			   CigarAlignment &altToRef,bool strict)
  : projected(projected), signals(signals), model(model), refSeq(refSeq), 
    refSeqStr(refSeqStr), altSeq(altSeq), altSeqStr(altSeqStr), G(NULL), 
    changes(false), altToRef(altToRef)
{
  if(!buildGraph(strict)) G=NULL;
}



ACEplus_Vertex *GraphBuilder::newVertex(const String &substrate,
					SignalType type,int begin,int end,
					double score,Strand strand,int ID,
					bool denovo)
{
  LightVertex *exists=G->vertexExists(substrate,strand,begin,end,type);
  if(exists) {
    if(denovo) dynamic_cast<ACEplus_Vertex*>(exists)->setDeNovo();
    return NULL;
  }
  SignalSensor *sensor=model.signalSensors->findSensor(type);
  const double rawScore=score;
  String signalStr;
  if(sensor) {
    int offset=sensor->getConsensusOffset();
    int contextLen=sensor->getContextWindowLength();
    int windowBegin=begin-offset;
    ContentSensor *bg=model.contentSensors->getSpliceBackground();
    double bgScore=bg->scoreSubsequence(altSeq,altSeqStr,windowBegin,
					contextLen,0);
    score-=bgScore; 
    signalStr=SignalPrinter::print(*sensor,windowBegin,altSeqStr);
  }
  ACEplus_Vertex *v=
    new ACEplus_Vertex(substrate,type,begin,end,score,strand,ID);
  v->setRawScore(rawScore);
  v->setSeq(signalStr);
  if(denovo) v->setDeNovo();
  return v;
}



ACEplus_Edge *GraphBuilder::newEdge(const String &substrate,ContentType type,
				 LightVertex *from,LightVertex *to,
				 int begin,int end,Strand strand,int ID)
{
  //cout<<"newEdge("<<type<<")"<<endl;
  if(end-begin<1) {
    //cout<<*G<<endl;
    cout<<"edge length < 1 : "<<begin<<"-"<<end<<" "<<substrate<<" "
	<<type<<endl;
    return NULL;
    //INTERNAL_ERROR
  }

  if(G->edgeExists(substrate,strand,begin,end,type))  {
    cout<<"edge exists -- ignoring: "<<begin<<"-"<<end<<" "<<type<<endl;
    return NULL;
  }
  //cout<<"returning new edge"<<endl;
  ACEplus_Edge *edge=
    new ACEplus_Edge(substrate,type,from,to,begin,end,strand,ID);
  edge->setScore(0);
  return edge;
}



LightGraph *GraphBuilder::getGraph()
{
  return G;
}



double GraphBuilder::scoreSignal(SignalType type,int pos,const Sequence &seq,
				const String &str)
{
  SignalSensor *sensor=model.signalSensors->findSensor(type);
  if(!sensor) return 0.0; // TSS and TES
  if(!sensor->consensusOccursAt(str,pos)) return NEGATIVE_INFINITY;
  int contextLen=sensor->getContextWindowLength();
  int offset=sensor->getConsensusOffset();
  int consensusLen=sensor->getConsensusLength();
  int begin=pos-offset;
  double score=sensor->getLogP(seq,str,begin);
  return score;
}



void GraphBuilder::inspectEdgeScore(ACEplus_Edge *edge,const String &label)
{
  //if(edge->getLength()<0) throw String("negative ")+edge->getLength()+"="+edge->getBegin()+"-"+edge->getEnd();
  const String id=projected.getGeneId();
  const PrefixSumArray &psa=model.contentSensors->getPSA(edge->getType());
  cout<<"EDGE\t"<<label<<"\t"<<id<<"\tL="<<edge->getLength()<<"\tLLR="
      <<psa.getInterval(edge->getBegin(),edge->getEnd())<<"\t";
  //for(int i=edge->getBegin() ; i<edge->getEnd() ; ++i) 
  //  cout<<psa.getInterval(i,i+1)<<"\t";
  cout<<endl;
}



double GraphBuilder::scoreEdge(ACEplus_Edge *edge)
{
  // Get emission probability
  if(edge->getLength()<1) {
    edge->setScore(NEGATIVE_INFINITY);
    return;
  }
  double score=model.contentSensors->score(edge->getType(),edge->getBegin(),
					  edge->getEnd());
  if(!isFinite(score)) {
    cerr<<*edge<<endl;
    cerr<<"edge emission score = "<<score<<" begin="<<edge->getBegin()<<
      " end="<<edge->getEnd()<<" type="<<edge->getType()<<endl;
    INTERNAL_ERROR; }
  const StructureChange change=edge->getChange();
  /*if(change.intronRetention || change.crypticExon)
    inspectEdgeScore(edge,change.intronRetention ? "intron-retention" :
    "denovo-exon");*/

  // Get duration probability
  const int L=edge->getLength();
  DiscreteDistribution *distr=NULL;
  if(edge->isExon()) distr=model.exonLengthDistr;
  else if(edge->isIntron()) distr=model.intronLengthDistr;
  else if(edge->isIntergenic()) distr=model.intergenicLengthDistr;
  else throw "Unknown edge type in GraphBuilder::scoreEdge()";
  const double durationScore=distr->getLogP(L);
  if(!isFinite(durationScore)) {
    cerr<<*edge<<endl;
    cerr<<"length="<<L<<", duration score="<<durationScore<<endl;
    INTERNAL_ERROR;
  }
  score+=durationScore;

  // Get transition probability
  const int numChoices=edge->getLeft()->getEdgesOut().size();
  //double uniform=1.0/numChoices;
  //const double transProb=log(uniform);
  SignalType fromType=edge->getLeft()->getSignalType();
  SignalType toType=edge->getRight()->getSignalType();
  const double transProb=model.transitions->getLogP(fromType,toType);
  if(!isFinite(transProb)) {
    cerr<<*edge<<endl;
    cerr<<"trans="<<transProb<<" from="<<fromType<<" to="<<toType<<endl;
    INTERNAL_ERROR;
  }
  score+=transProb; //###

  if(numChoices==0) {
    cout<<*edge->getLeft()<<endl;
    cout<<*edge<<endl;
    cerr<<"Zero edges out!"<<endl;
    INTERNAL_ERROR; }

  // Apply hard filter to reduce false positives
  if(change.intronRetention && score<model.minIntronRetentionLLR)
    score=NEGATIVE_INFINITY;
  else if(change.crypticExon && score<model.minDeNovoExonLLR)
    score=NEGATIVE_INFINITY;

  // Store score in edge
  edge->setScore(score);
}



void GraphBuilder::getContextWindow(LightVertex *v,int &begin,int &end)
{
  SignalSensor *sensor=model.signalSensors->findSensor(v->getType());
  if(!sensor) { begin=v->getBegin(); end=v->getEnd(); return; }
  int consensusPos=v->getBegin();
  begin=consensusPos-sensor->getConsensusOffset();
  end=begin+sensor->getContextWindowLength();
}



void GraphBuilder::getSignalWindow(SignalType type,int consensusPos,
				   int &windowBegin,int &windowEnd,
				   int &consensusLen)
{
  if(type==TSS || type==TES) {
    windowBegin=windowEnd=consensusPos;
    consensusLen=0;
    return; }
  SignalSensor *sensor=model.signalSensors->findSensor(type);
  if(!sensor) INTERNAL_ERROR;
  int contextLen=sensor->getContextWindowLength();
  int offset=sensor->getConsensusOffset();
  consensusLen=sensor->getConsensusLength();
  windowBegin=consensusPos-offset;
  windowEnd=windowBegin+contextLen;
}



ContentType GraphBuilder::getContentType(SignalType from,SignalType to)
{
  switch(from) {
  case LEFT_TERMINUS: return INTERGENIC;
  case TSS:           return EXON;
  case GT:            return INTRON;
  case AG:            return EXON;
  case TES:           return INTERGENIC;
  }
  INTERNAL_ERROR;
}



bool GraphBuilder::buildTranscriptGraph()
{
  Strand strand=projected.getStrand();
  if(strand!=FORWARD_STRAND) INTERNAL_ERROR;
  const String &substrate=projected.getSubstrate();
  G=new LightGraph(substrate,altSeq.getLength());
  int numSignals=signals.numSignals();

  // Create vertices
  ACEplus_Vertex *v=newVertex(substrate,LEFT_TERMINUS,0,0,0.0,strand,0);
  v->setAnnotated(true); v->setBroken(false);
  ACEplus_Vertex *leftTerminus=v;
  G->addVertex(v);
  for(int i=0 ; i<numSignals ; ++i) {
    const TranscriptSignal &signal=signals[i];
    SignalType type=signal.getType();
    int pos=signal.getPos();
    bool broken=signal.isBroken();
    double score=scoreSignal(type,pos,altSeq,altSeqStr);
    int begin, end, consensusLen;
    getSignalWindow(type,pos,begin,end,consensusLen);
    ACEplus_Vertex *v=newVertex(substrate,type,pos,pos+consensusLen,score,
				   strand,i+1);
#ifdef SANITY_CHECKS
    if(!v) INTERNAL_ERROR;
#endif
    v->setBroken(signal.isBroken());
    v->setAnnotated(true);
    G->addVertex(v);
  }
  int L=altSeq.getLength();
  v=newVertex(substrate,RIGHT_TERMINUS,L,L,0.0,strand,numSignals+1);
#ifdef SANITY_CHECKS
  if(!v) INTERNAL_ERROR;
#endif
  v->setAnnotated(true); v->setBroken(false);
  ACEplus_Vertex *rightTerminus=v;
  G->addVertex(v);

  // Create edges
  const int numVertices=G->getNumVertices();
  for(int i=0 ; i+1<numVertices ; ++i) {
    LightVertex *prev=G->getVertex(i), *next=G->getVertex(i+1);
    int begin, end, dummy;
    getContextWindow(prev,dummy,begin); getContextWindow(next,end,dummy);
    ContentType type=getContentType(prev->getType(),next->getType());
#ifdef SANITY_CHECKS
    if(prev->getType()==next->getType()) INTERNAL_ERROR;
#endif
    LightEdge *edge=newEdge(substrate,type,prev,next,begin,end,strand,i);
    if(!edge) {
      //ACEplus_Edge *edge=
      //GraphBuilder::linkVertices(leftTerminus,rightTerminus);
      return false;
      //INTERNAL_ERROR;
    }
    edge->setBroken(false);
    edge->setAnnotated(true);
    prev->addEdgeOut(edge); next->addEdgeIn(edge);
    G->addEdge(edge);
  }  
  return true;
}



void GraphBuilder::handleBrokenSites()
{
  // For each broken site, scan for cryptic sites nearby
  int numVertices=G->getNumVertices();
  for(int i=0 ; i<numVertices ; ++i) {
    LightVertex *v=G->getVertex(i);
    if(v->isBroken()) handleBrokenSite(v);
  }

  // Delete broken sites
  Set<LightVertex*> verticesToDelete;
  Set<LightEdge*> edgesToDelete;
  for(int i=0 ; i<numVertices ; ++i) {
    LightVertex *v=G->getVertex(i);
    if(!v->isBroken()) continue;
    verticesToDelete+=v;
    Vector<LightEdge*> &in=v->getEdgesIn(), &out=v->getEdgesOut();
    for(Vector<LightEdge*>::iterator cur=in.begin(), end=in.end() ; 
	cur!=end ; ++cur) edgesToDelete+=*cur;
    for(Vector<LightEdge*>::iterator cur=out.begin(), end=out.end() ; 
	cur!=end ; ++cur) edgesToDelete+=*cur;
    G->dropVertex(i);
  }
  for(Set<LightEdge*>::iterator cur=edgesToDelete.begin(), end=
	edgesToDelete.end() ; cur!=end ; ++cur) {
    LightEdge *e=*cur;
    e->getLeft()->dropEdgeOut(e); e->getRight()->dropEdgeIn(e);
    G->dropEdge(e->getID());
    delete e;
  }
  for(Set<LightVertex*>::iterator cur=verticesToDelete.begin(), end=
	verticesToDelete.end() ; cur!=end ; ++cur) delete *cur;
  G->deleteNullVertices(); G->deleteNullEdges();
}



void GraphBuilder::handleBrokenSite(LightVertex *v)
{
  if(model.allowExonSkipping) handleBrokenSite_skipping(v);
  if(model.allowIntronRetention) handleBrokenSite_retention(v);
  if(model.allowCrypticSites) handleBrokenSite_cryptic(v);
}




void GraphBuilder::handleBrokenSite_cryptic(LightVertex *v)
{
  // Prepare to do the scan
  SignalType type=v->getType();
  const String &substrate=v->getSubstrate();
  const Strand strand=v->getStrand();
  SignalSensor *sensor=model.signalSensors->findSensor(type);
  const double cutoff=sensor->getCutoff();
  const int contextWindowLen=sensor->getContextWindowLength();
  const int consensusOffset=sensor->getConsensusOffset();
  const int consensusLen=sensor->getConsensusLength();
  int scanBegin=v->getBegin()-model.MAX_SPLICE_SHIFT-consensusOffset;
  int scanEnd=v->getEnd()+model.MAX_SPLICE_SHIFT-consensusLen-consensusOffset;
  if(scanBegin<0) scanBegin=0;
  if(scanEnd>altSeq.getLength()-contextWindowLen)
    scanEnd=altSeq.getLength()-contextWindowLen;
#ifdef SANITY_CHECKS
  if(scanEnd<0) INTERNAL_ERROR;
#endif
  Vector<LightEdge*> &in=v->getEdgesIn(), &out=v->getEdgesOut();
  int nextVertexID=G->getNumVertices(), nextEdgeID=G->getNumEdges();

  // Scan in a window around the broken site
  for(int pos=scanBegin ; pos<scanEnd ; ++pos) {
    const int consensusPos=pos+consensusOffset;
    if(consensusPos==v->getBegin()) continue;
    if(sensor->consensusOccursAt(altSeqStr,consensusPos)) {
      double score=sensor->getLogP(altSeq,altSeqStr,pos);
      if(score<cutoff) continue;
      ACEplus_Vertex *v=newVertex(substrate,type,consensusPos,
			       consensusPos+consensusLen,score,
			       strand,nextVertexID++);
      if(!v) continue;
      v->setBroken(false);  v->setAnnotated(false);
      v->setThreshold(sensor->getCutoff());
      v->setRefScore(0); // ###
      String signalStr=SignalPrinter::print(*sensor,pos,altSeqStr);
      v->setSeq(signalStr);
      G->addVertex(v);
      for(Vector<LightEdge*>::iterator cur=in.begin(), end=in.end() ; 
	  cur!=end ; ++cur) {
	LightEdge *edge=*cur; LightVertex *left=edge->getLeft();
	int begin, end, dummy;
	getContextWindow(left,dummy,begin); getContextWindow(v,end,dummy);
#ifdef SANITY_CHECKS
	if(left->getType()==v->getType()) INTERNAL_ERROR;
#endif
	ACEplus_Edge *new_Edge=newEdge(substrate,edge->getType(),left,v,
				   begin,end,strand,nextEdgeID++);
	if(!new_Edge) continue;
	new_Edge->getChange().crypticSite=true;
	new_Edge->setBroken(false);
	left->addEdgeOut(new_Edge);
	v->addEdgeIn(new_Edge); G->addEdge(new_Edge);
      }
      for(Vector<LightEdge*>::iterator cur=out.begin(), end=out.end() ; 
	  cur!=end ; ++cur) {
	LightEdge *edge=*cur; LightVertex *right=edge->getRight();
	int begin, end, dummy;
	getContextWindow(v,dummy,begin); getContextWindow(right,end,dummy);
#ifdef SANITY_CHECKS
	if(v->getType()==right->getType()) INTERNAL_ERROR;
#endif
	ACEplus_Edge *new_Edge=newEdge(substrate,edge->getType(),v,right,
				   begin,end,strand,nextEdgeID++);
	if(!new_Edge) continue;
	new_Edge->getChange().crypticSite=true;
	new_Edge->setBroken(false);
	right->addEdgeIn(new_Edge);
	v->addEdgeOut(new_Edge); G->addEdge(new_Edge);
      }
    }
  }
}



LightEdge *GraphBuilder::getAnnotatedEdge(Vector<LightEdge*> &edges)
{
  for(Vector<LightEdge*>::iterator cur=edges.begin(), end=edges.end() ;
      cur!=end ; ++cur) {
    LightEdge *edge=*cur;
    if(edge->isAnnotated()) return edge; }
  INTERNAL_ERROR;
}



void GraphBuilder::handleBrokenSite_skipping(LightVertex *v)
{
  Vector<LightEdge*> &in=v->getEdgesIn(), &out=v->getEdgesOut();
  if(in.size()==0 || out.size()==0) return;
  //if(in[0]->getType()==EXON) { handleSkippingLeft(v); return; }
  //if(out[0]->getType()==EXON) { handleSkippingRight(v); return; }
  if(getAnnotatedEdge(in)->getType()==EXON) handleSkippingLeft(v);
  else if(getAnnotatedEdge(out)->getType()==EXON) handleSkippingRight(v);
}



void GraphBuilder::handleSkippingLeft(LightVertex *v)
{
  Vector<LightEdge*> &in=v->getEdgesIn(), &out=v->getEdgesOut();
  //LightEdge *leftEdge=in[0], *rightEdge=out[0];
  LightEdge *leftEdge=getAnnotatedEdge(in), *rightEdge=getAnnotatedEdge(out);
#ifdef SANITY_CHECKS
  if(!leftEdge || !rightEdge) INTERNAL_ERROR;
#endif
  const String &substrate=leftEdge->getSubstrate();
  Strand strand=leftEdge->getStrand();
  if(rightEdge->getType()!=INTRON) return;
  LightVertex *to=rightEdge->getRight();
  Vector<LightEdge*> &leftIntrons=leftEdge->getLeft()->getEdgesIn();
  if(leftIntrons.size()==0 || leftIntrons[0]->getType()!=INTRON) return;
  for(Vector<LightEdge*>::iterator cur=leftIntrons.begin(), 
	end=leftIntrons.end() ; cur!=end ; ++cur) {
    LightVertex *from=(*cur)->getLeft();
    int begin, end, dummy;
    getContextWindow(from,dummy,begin); getContextWindow(to,end,dummy);
#ifdef SANITY_CHECKS
    if(from->getType()==to->getType()) INTERNAL_ERROR;
#endif
    ACEplus_Edge *new_Edge=
      newEdge(substrate,INTRON,from,to,begin,end,strand,
	      G->getNumEdges());
    if(!new_Edge) continue;
    new_Edge->getChange().exonSkipping=true;
    new_Edge->setBroken(false);
    from->addEdgeOut(new_Edge); to->addEdgeIn(new_Edge);
    G->addEdge(new_Edge);
  }
}



void GraphBuilder::handleSkippingRight(LightVertex *v)
{
  Vector<LightEdge*> &in=v->getEdgesIn(), &out=v->getEdgesOut();
  LightEdge *leftEdge=getAnnotatedEdge(in), *rightEdge=getAnnotatedEdge(out);
#ifdef SANITY_CHECKS
  if(!leftEdge || !rightEdge) INTERNAL_ERROR;
#endif
  const String &substrate=leftEdge->getSubstrate();
  Strand strand=leftEdge->getStrand();
  if(leftEdge->getType()!=INTRON) return;
  LightVertex *from=leftEdge->getLeft();
  Vector<LightEdge*> &rightIntrons=rightEdge->getRight()->getEdgesOut();
  if(rightIntrons.size()==0 || rightIntrons[0]->getType()!=INTRON) return;
  for(Vector<LightEdge*>::iterator cur=rightIntrons.begin(), 
	end=rightIntrons.end() ; cur!=end ; ++cur) {
    LightVertex *to=(*cur)->getRight();
    int begin, end, dummy;
    getContextWindow(from,dummy,begin); getContextWindow(to,end,dummy);
#ifdef SANITY_CHECKS
    if(from->getType()==to->getType()) INTERNAL_ERROR;
#endif
    ACEplus_Edge *new_Edge=
      newEdge(substrate,INTRON,from,to,begin,end,strand,
	      G->getNumEdges());
    if(!new_Edge) continue;
    new_Edge->getChange().exonSkipping=true;
    new_Edge->setBroken(false);
    from->addEdgeOut(new_Edge); to->addEdgeIn(new_Edge);
    G->addEdge(new_Edge);
  }
}



void GraphBuilder::handleBrokenSite_retention(LightVertex *v)
{
  Vector<LightEdge*> &in=v->getEdgesIn(), &out=v->getEdgesOut();
  if(in.size()==0 || out.size()==0) return;
  if(getAnnotatedEdge(in)->getType()==INTRON) handleRetentionLeft(v);
  else if(getAnnotatedEdge(out)->getType()==INTRON) handleRetentionRight(v);
}



void GraphBuilder::handleRetentionLeft(LightVertex *v)
{
  Vector<LightEdge*> &in=v->getEdgesIn(), &out=v->getEdgesOut();
  LightEdge *leftEdge=getAnnotatedEdge(in), *rightEdge=getAnnotatedEdge(out);
#ifdef SANITY_CHECKS
  if(!leftEdge || !rightEdge) INTERNAL_ERROR;
#endif
  const String &substrate=leftEdge->getSubstrate();
  Strand strand=leftEdge->getStrand();
  LightVertex *to=rightEdge->getRight();
  Vector<LightEdge*> &leftExons=leftEdge->getLeft()->getEdgesIn();
  for(Vector<LightEdge*>::iterator cur=leftExons.begin(), 
	end=leftExons.end() ; cur!=end ; ++cur) {
    LightVertex *from=(*cur)->getLeft();
    int begin, end, dummy;
    getContextWindow(from,dummy,begin); getContextWindow(to,end,dummy);
#ifdef SANITY_CHECKS
    if(from->getType()==to->getType()) INTERNAL_ERROR;
#endif
    if(end-begin>model.maxIntronRetentionLen) continue;
    ACEplus_Edge *new_Edge=
      newEdge(substrate,EXON,from,to,begin,end,strand,G->getNumEdges());
    if(!new_Edge) continue;
    new_Edge->getChange().intronRetention=true;
    from->addEdgeOut(new_Edge); to->addEdgeIn(new_Edge);
    new_Edge->setBroken(false);
    G->addEdge(new_Edge);
  }
}



void GraphBuilder::handleRetentionRight(LightVertex *v)
{
  Vector<LightEdge*> &in=v->getEdgesIn(), &out=v->getEdgesOut();
  LightEdge *leftEdge=getAnnotatedEdge(in), *rightEdge=getAnnotatedEdge(out);
#ifdef SANITY_CHECKS
  if(!leftEdge || !rightEdge) INTERNAL_ERROR;
#endif
  const String &substrate=leftEdge->getSubstrate();
  Strand strand=leftEdge->getStrand();
  LightVertex *from=leftEdge->getLeft();
  Vector<LightEdge*> &rightExons=rightEdge->getRight()->getEdgesOut();
  for(Vector<LightEdge*>::iterator cur=rightExons.begin(), 
	end=rightExons.end() ; cur!=end ; ++cur) {
    LightVertex *to=(*cur)->getRight();
    int begin, end, dummy;
    getContextWindow(from,dummy,begin); getContextWindow(to,end,dummy);
#ifdef SANITY_CHECKS
    if(from->getType()==to->getType()) INTERNAL_ERROR;
#endif
    if(end-begin>model.maxIntronRetentionLen) continue;
    ACEplus_Edge *new_Edge=
      newEdge(substrate,EXON,from,to,begin,end,strand,G->getNumEdges());
    if(!new_Edge) continue;
    new_Edge->getChange().intronRetention=true;
    from->addEdgeOut(new_Edge); to->addEdgeIn(new_Edge);
    new_Edge->setBroken(false);
    G->addEdge(new_Edge);
  }
}



void GraphBuilder::leftSweep(LightGraph &G,
			     Array1D<int> &vertexCounts,
			     Array1D<int> &edgeCounts)
{
  if(G.getNumVertices()==0) throw "empty graph";
  LightVertex *left=G.getVertex(0);
  Stack<LightVertex*> todo;
  todo.push(left);
  while(!todo.isEmpty()) {
    LightVertex *v=todo.pop();
    ++vertexCounts[v->getID()];
    Vector<LightEdge*> &edges=v->getEdgesOut();
    for(Vector<LightEdge*>::iterator cur=edges.begin(), end=edges.end() ;
	cur!=end ; ++cur) {
      LightEdge *edge=*cur;
      ++edgeCounts[edge->getID()];
      if(edge->getRight()->getBegin()<v->getBegin()) {
	cout<<*v<<"\n"<<*edge->getRight()<<endl;
	INTERNAL_ERROR;
      }
      if(vertexCounts[edge->getRight()->getID()]==0)
	todo.push(edge->getRight());
    }
  }
}



void GraphBuilder::rightSweep(LightGraph &G,
			      Array1D<int> &vertexCounts,
			      Array1D<int> &edgeCounts)
{
  const int numVertices=G.getNumVertices();
  if(numVertices==0) throw "empty graph";
  LightVertex *right=G.getVertex(numVertices-1);
  Stack<LightVertex*> todo;
  todo.push(right);
  while(!todo.isEmpty()) {
    LightVertex *v=todo.pop();
    ++vertexCounts[v->getID()];
    Vector<LightEdge*> &edges=v->getEdgesIn();
    for(Vector<LightEdge*>::iterator cur=edges.begin(), end=edges.end() ;
	cur!=end ; ++cur) {
      LightEdge *edge=*cur;
      ++edgeCounts[edge->getID()];
      if(vertexCounts[edge->getLeft()->getID()]==0)
	todo.push(edge->getLeft());
    }
  }
}



void GraphBuilder::deleteUnreachable(LightGraph &G,
				     Array1D<int> &vertexLeftCounts,
				     Array1D<int> &vertexRightCounts,
				     Array1D<int> &edgeLeftCounts,
				     Array1D<int> &edgeRightCounts)
{
  for(int i=0 ; i<vertexLeftCounts.size() ; ++i) 
    if(vertexLeftCounts[i]==0 || vertexRightCounts[i]==0) G.dropVertex(i);
  for(int i=0 ; i<edgeLeftCounts.size() ; ++i) 
    if(edgeLeftCounts[i]==0 || edgeRightCounts[i]==0) {
      LightEdge *edge=G.getEdge(i);
      edge->getLeft()->dropEdgeOut(edge);
      edge->getRight()->dropEdgeIn(edge);
      G.dropEdge(i);
    }
}



void GraphBuilder::pruneUnreachable(LightGraph &G)
{
  G.sort();
  Array1D<int> vertexLeftCounts(G.getNumVertices());
  Array1D<int> edgeLeftCounts(G.getNumEdges());
  Array1D<int> vertexRightCounts(G.getNumVertices());
  Array1D<int> edgeRightCounts(G.getNumEdges());
  vertexLeftCounts.setAllTo(0);
  vertexRightCounts.setAllTo(0);
  edgeLeftCounts.setAllTo(0);
  edgeRightCounts.setAllTo(0);
  //cout<<"leftSweep"<<endl;
  leftSweep(G,vertexLeftCounts,edgeLeftCounts);
  //cout<<"rightSweep"<<endl;
  rightSweep(G,vertexRightCounts,edgeRightCounts);
  //cout<<"deleteUnreachable"<<endl;
  deleteUnreachable(G,vertexLeftCounts,vertexRightCounts,
		    edgeLeftCounts,edgeRightCounts);
  //cout<<"delete null vertices"<<endl;
  G.deleteNullVertices();
  //cout<<"delete null edge"<<endl;
  G.deleteNullEdges();
  //cout<<"sort"<<endl;
  G.sort();
}



bool GraphBuilder::mapped() const
{
  return !changes;
}



void GraphBuilder::findVariants(Vector<Interval> &variants)
{
  Vector<int> positions;
  const int L=altSeq.getLength();
  for(int pos=0 ; pos<L ; ++pos) {
    const int to=altToRef[pos];
    if(to==CIGAR_UNDEFINED ||              // Insertion
       pos+1<L && altToRef[pos+1]!=to+1 || // Deletion
       altSeqStr[pos]!=refSeqStr[to])      // Substitution
      positions.push_back(pos); }
  const int N=positions.size();
  int begin=-1, end=-1;
  for(int i=0 ; i<N ; ++i) {
    const int pos=positions[i];
    if(begin<0) { begin=pos; end=pos+1; }
    else if(pos==end) ++end;
    else { 
      variants.push_back(Interval(begin,end));
      begin=end=-1; }
  }
  if(begin>=0) variants.push_back(Interval(begin,end));
}



void GraphBuilder::handleDeNovoSites()
{
  cout<<"handleDeNovoSites()"<<endl;

  // First, find all variants
  if(variants.isEmpty()) findVariants(variants);
  cout<<variants.size()<<" variants found"<<endl;

  // Now scan for de novo sites overlapping any variant
  Vector<ACEplus_Vertex*> newVertices;
  scanDeNovo(*model.signalSensors->donorSensor,variants,newVertices);
  cout<<newVertices.size()<<" new vertices after scanning for donors"<<endl;
  scanDeNovo(*model.signalSensors->acceptorSensor,variants,newVertices);
  cout<<newVertices.size()<<" after scanning for acceptors"<<endl;

  // Now link the new vertices to other vertices
  cout<<G->getNumVertices()<<" vertices before linking"<<endl;
  linkDeNovoVertices(newVertices);

  // If cryptic exons are enabled, scan for partner signals to complete
  // new exons
  cout<<G->getNumVertices()<<" vertices before adding cryptic exons"<<endl;
  if(model.allowCrypticExons) addCrypticExons(newVertices);
  cout<<newVertices.size()<<" after scanning for cryptic exons"<<endl;
  cout<<G->getNumVertices()<<" vertices after adding cryptic exons"<<endl;
}



void dumpIntervals(const Vector<Interval> &windows)
{
  for(Vector<Interval>::const_iterator cur=windows.begin(), end=windows.end() ;
      cur!=end ; ++cur)
    cout<<*cur<<endl;
}



void GraphBuilder::scanDeNovo(SignalSensor &sensor,Vector<Interval> &variants,
			      Vector<ACEplus_Vertex*> &newVertices)
{
  // PRECONDITION: variants are sorted by increasing absolute position

  // First, get a set of windows to run the sensor on
  Vector<Interval> windows;
  getVariantWindows(windows,sensor,variants);
  //cout<<windows.size()<<" variant windows before coalescing:"<<endl;
  //dumpIntervals(windows);

  // Coalesce windows, in case any overlap
  coalesceWindows(windows);
  //cout<<windows.size()<<" variant windows after coalescing:"<<endl;
  //dumpIntervals(windows);

  // Perform signal sensing in all windows
  deNovoSignalSensing(sensor,windows,newVertices);
}



void GraphBuilder::getVariantWindows(Vector<Interval> &windows,
				     SignalSensor &sensor,
				     Vector<Interval> &variants)
{
  const int sensorLen=sensor.getContextWindowLength();
  for(Vector<Interval>::iterator cur=variants.begin(), end=variants.end() ;
      cur!=end ; ++cur) {
    const Interval &variant=*cur;
    int begin=variant.getBegin()-sensorLen+1;
    if(begin<0) begin=0;
    int end=variant.getEnd();
    windows.push_back(Interval(begin,end));
  }
}



void GraphBuilder::coalesceWindows(Vector<Interval> &windows)
{
  int numWindows=windows.size();
  for(int i=0 ; i+1<numWindows ; ++i) {
    if(windows[i].overlaps(windows[i+1])) {
      windows[i].setEnd(windows[i+1].getEnd());
      windows.cut(i+1);
      --i; --numWindows;
    }
  }
}



void GraphBuilder::deNovoSignalSensing(SignalSensor &sensor,
				       Vector<Interval> &windows,
				       Vector<ACEplus_Vertex*> &newVertices)
{
  const double threshold=sensor.getCutoff();
  const int consensusOffset=sensor.getConsensusOffset();
  const int consensusLen=sensor.getConsensusLength();
  const String &substrate=projected.getSubstrate();
  const SignalType signalType=sensor.getSignalType();
  const Strand strand=projected.getStrand();
  const int numWindows=windows.size();
  const int sensorLen=sensor.getContextWindowLength();
  for(int i=0 ; i<numWindows ; ++i) {
    const Interval &scanWindow=windows[i];
    for(int pos=scanWindow.getBegin() ; pos<scanWindow.getEnd() ; ++pos) {
      if(pos+sensorLen>altSeq.getLength()) continue;
      if(sensor.consensusOccursAt(altSeqStr,pos+consensusOffset)) {
	const double altScore=sensor.getLogP(altSeq,altSeqStr,pos);
	if(altScore<threshold) continue;
	const int refPos=altToRef[pos];
	if(refPos!=CIGAR_UNDEFINED && refPos+sensorLen<=refSeq.getLength() &&
	   sensor.consensusOccursAt(refSeqStr,refPos+consensusOffset)) {
	  const double refScore=sensor.getLogP(refSeq,refSeqStr,refPos);
	  if(refScore>=threshold) continue; // ref already had a signal there
	  if(altScore-refScore<log(2)) continue; // less than a 2-fold increase
	}
	if(altScore<threshold+log(2)) continue; // must be >= 2*threshold
	ACEplus_Vertex *v=newVertex(substrate,signalType,pos+consensusOffset,
				    pos+consensusOffset+consensusLen,
				    altScore,strand,G->getNumVertices(),
				    true);
	if(!v) continue;
	G->addVertex(v);
	newVertices.push_back(v);
	v->setThreshold(threshold);
      }
    }
  }
}



void GraphBuilder::linkDeNovoVertices(Vector<ACEplus_Vertex*> &newVertices)
{
  G->sort();
  LightGraph &graph=*G;
  for(Vector<ACEplus_Vertex*>::iterator cur=newVertices.begin(), 
	end=newVertices.end() ; cur!=end ; ++cur) {
    ACEplus_Vertex *v=*cur;
    const int id=v->getID();
#ifdef SANITY_CHECKS
    if(graph.getVertex(id)!=v) INTERNAL_ERROR;
#endif
    linkDeNovoLeft(v,id);
    linkDeNovoRight(v,id);
  }
}



void GraphBuilder::getLinkTypes(SignalType t,Direction dir,
				Set<SignalType> &types)
{
  switch(dir) {
  case DIR_LEFT:
    switch(t) {
    case GT: types.insert(TSS); types.insert(AG); break;
    case AG: types.insert(GT); break; }
    break;
  case DIR_RIGHT:
    switch(t) {
    case GT: types.insert(AG); break;
    case AG: types.insert(GT); types.insert(TES); break; }
    break;
  default: INTERNAL_ERROR;
  }
}



void GraphBuilder::linkDeNovoLeft(ACEplus_Vertex *v,int id)
{
  const String &substrate=projected.getSubstrate();
  Strand strand=projected.getStrand();
  LightGraph &graph=*G;
  Set<SignalType> linkTypes;
  getLinkTypes(v->getType(),DIR_LEFT,linkTypes);
  for(int i=id-1 ; i>=0 ; --i) {
    LightVertex *w=graph.getVertex(i);
    if(linkTypes.isMember(w->getType())) {
      const int edgeID=graph.getNumEdges();
      int begin, end, dummy;
      getContextWindow(w,dummy,begin); getContextWindow(v,end,dummy);
      ContentType type=getContentType(w->getType(),v->getType());
#ifdef SANITY_CHECKS
      if(w->getType()==v->getType()) INTERNAL_ERROR;
#endif
      //if(::isExon(type) && v->distanceTo(*w)>model.maxDeNovoExonLen) continue;
      ACEplus_Edge *edge=
	newEdge(substrate,type,w,v,begin,end,strand,edgeID);
      if(!edge) continue;
      edge->getChange().deNovoSite=true;
      edge->setBroken(false);
      graph.addEdge(edge);
      w->addEdgeOut(edge); v->addEdgeIn(edge);
      if(w->isAnnotated()) break;
    }
  }
}



void GraphBuilder::linkDeNovoRight(ACEplus_Vertex *v,int id)
{
  const String &substrate=projected.getSubstrate();
  Strand strand=projected.getStrand();
  LightGraph &graph=*G;
  Set<SignalType> linkTypes;
  getLinkTypes(v->getType(),DIR_RIGHT,linkTypes);
  const int numVertices=graph.getNumVertices();
  for(int i=id+1 ; i<numVertices ; ++i) {
    LightVertex *w=graph.getVertex(i);
    if(linkTypes.isMember(w->getType())) {
      const int edgeID=graph.getNumEdges();
      int begin, end, dummy;
      getContextWindow(v,dummy,begin); getContextWindow(w,end,dummy);
      ContentType type=getContentType(v->getType(),w->getType());
#ifdef SANITY_CHECKS
      if(v->getType()==w->getType()) INTERNAL_ERROR;
#endif
      //if(::isExon(type) && v->distanceTo(*w)>model.maxDeNovoExonLen) continue;
      ACEplus_Edge *edge=
	newEdge(substrate,type,v,w,begin,end,strand,edgeID);
      if(!edge) continue;
      edge->getChange().deNovoSite=true;
      edge->setBroken(false);
      graph.addEdge(edge);
      v->addEdgeOut(edge); w->addEdgeIn(edge);
      if(w->isAnnotated()) break;
    }
  }
}



bool GraphBuilder::allVerticesAreAnnotated()
{  
  LightGraph &graph=*G;
  const int numVertices=graph.getNumVertices();
  //cout<<numVertices<<" vertices"<<endl;
  for(int i=0 ; i<numVertices ; ++i) {
    LightVertex *v=graph.getVertex(i);
    if(!v->isAnnotated()) return false;
  }
  //cout<<"returning true"<<endl;
  return true;
}



void GraphBuilder::addCrypticExons(Vector<ACEplus_Vertex*> &newVertices)
{
  G->sort();
  cout<<newVertices.size()<<" new vertices"<<endl;
  for(Vector<ACEplus_Vertex*>::iterator cur=newVertices.begin(), end=
	newVertices.end() ; cur!=end ; ++cur) {
    ACEplus_Vertex *v=*cur;
    const SignalType t=v->getType();
    const int id=v->getID();
    //cout<<"scanning..."<<t<<endl;
    if(isDonor(t)) scanCrypticExonLeft(id);
    else if(isAcceptor(t)) scanCrypticExonRight(id);
    else INTERNAL_ERROR;
    //cout<<"scan finished"<<endl;
  }
}



LightVertex *GraphBuilder::findAnnotatedVertexLeftOf(int of)
{
  for(int i=of-1 ; i>=0 ; --i) {
    LightVertex *v=G->getVertex(i);
    if(v->isAnnotated()) return v;
  }
  return NULL;
}



LightVertex *GraphBuilder::findAnnotatedVertexRightOf(int of)
{
  const int numVertices=G->getNumVertices();
  for(int i=of+1 ; i<numVertices ; ++i) {
    LightVertex *v=G->getVertex(i);
    if(v->isAnnotated()) return v;
  }
  return NULL;
}



void GraphBuilder::findMateSignals(SignalSensor &sensor,
				   const Interval &scanWindow,
				   Vector<ACEplus_Vertex*> &newVertices)
{
  /* Once a variant creates a new splice site in an intron, this function
     finds a mate splice site that can complete a cryptic exon */
  const int consensusOffset=sensor.getConsensusOffset();
  const int consensusLen=sensor.getConsensusLength();
  const String &substrate=projected.getSubstrate();
  const SignalType signalType=sensor.getSignalType();
  const Strand strand=projected.getStrand();
  const int sensorLen=sensor.getContextWindowLength();
  for(int pos=scanWindow.getBegin() ; pos<scanWindow.getEnd()-sensorLen+1 ; 
      ++pos) {
    if(pos+sensorLen>altSeq.getLength()) continue;
    if(sensor.consensusOccursAt(altSeqStr,pos+consensusOffset)) {
      const double altScore=sensor.getLogP(altSeq,altSeqStr,pos);
      if(altScore<sensor.getCutoff()) continue;
      ACEplus_Vertex *v=newVertex(substrate,signalType,pos+consensusOffset,
			       pos+consensusOffset+consensusLen,
			       altScore,strand,G->getNumVertices());
      if(!v) continue;
      G->addVertex(v);
      newVertices.push_back(v);
      v->setThreshold(sensor.getCutoff());
    }
  }
}



ACEplus_Edge *GraphBuilder::linkVertices(LightVertex *left,LightVertex *right)
{
#ifdef SANITY_CHECKS
  if(left->getEnd()>=right->getBegin()) {
    cout<<*left<<"\n"<<*right<<endl;
    INTERNAL_ERROR;
  }
#endif
  int leftBegin, leftEnd, rightBegin, rightEnd;
  getContextWindow(left,leftBegin,leftEnd);
  getContextWindow(right,rightBegin,rightEnd);

  const String &substrate=projected.getSubstrate();
  const Strand strand=projected.getStrand();
  ContentType type=getContentType(left->getType(),right->getType());
  int edgeID=G->getNumEdges();
#ifdef SANITY_CHECKS
  if(left->getType()==right->getType()) INTERNAL_ERROR;
#endif
  ACEplus_Edge *edge=
    newEdge(substrate,type,left,right,leftEnd,rightBegin,strand,edgeID);
  if(!edge) return NULL;
  edge->setBroken(false);
  G->addEdge(edge);
  left->addEdgeOut(edge);
  right->addEdgeIn(edge);
  return edge;
}



/* This function scans left of a new donor site to search for a matching
   acceptor site that could complete a cryptic exon deep in an intron.
 */
void GraphBuilder::scanCrypticExonLeft(int of)
{
  // First, identify the scanning region
  int ofBegin, ofEnd, leftBegin, leftEnd;
  LightVertex *ofVertex=G->getVertex(of);
#ifdef SANITY_CHECKS
  if(ofVertex->getType()!=GT) INTERNAL_ERROR;
#endif
  getContextWindow(ofVertex,ofBegin,ofEnd);
  LightVertex *left=findAnnotatedVertexLeftOf(of);
#ifdef SANITY_CHECKS
  if(!left) INTERNAL_ERROR;
#endif
  if(left->getType()!=GT) return;
  getContextWindow(left,leftBegin,leftEnd);

  // Scan
  SignalSensor *sensor=model.signalSensors->findSensor(AG);

  Vector<ACEplus_Vertex*> newVertices;
  if(ofBegin-leftEnd>model.maxDeNovoExonLen)
    leftEnd=ofBegin-model.maxDeNovoExonLen;
  findMateSignals(*sensor,Interval(leftEnd,ofBegin),newVertices);
  
  // Link the new vertices into the graph
  for(Vector<ACEplus_Vertex*>::iterator cur=newVertices.begin(), end=
	newVertices.end() ; cur!=end ; ++cur) {
    LightVertex *newVertex=*cur;
    //int newBegin, newEnd;
    //getContextWindow(newVertex,newBegin,newEnd);

    // Link this new vertex to the de novo vertex, to create a cryptic exon
#ifdef SANITY_CHECKS
    if(newVertex->getType()==ofVertex->getType()) INTERNAL_ERROR;
#endif
    ACEplus_Edge *edge=linkVertices(newVertex,ofVertex);
    if(edge) edge->getChange().crypticExon=true;

    // Link this new vertex left to the nearest annotated vertex
#ifdef SANITY_CHECKS
    if(left->getType()==newVertex->getType()) INTERNAL_ERROR;
#endif
    linkVertices(left,newVertex);
  }
}



/* This function scans right of a new acceptor site to search for a matching
   donor site that could complete a cryptic exon deep in an intron.
 */
void GraphBuilder::scanCrypticExonRight(int of)
{
  // First, identify the scanning region
  int ofBegin, ofEnd, rightBegin, rightEnd;
  LightVertex *ofVertex=G->getVertex(of);
#ifdef SANITY_CHECKS
  if(ofVertex->getType()!=AG) INTERNAL_ERROR;
#endif
  getContextWindow(ofVertex,ofBegin,ofEnd);
  LightVertex *right=findAnnotatedVertexRightOf(of);
#ifdef SANITY_CHECKS
  if(!right) INTERNAL_ERROR;
#endif
  if(right->getType()!=AG) return;
  getContextWindow(right,rightBegin,rightEnd);

  // Scan
  SignalSensor *sensor=model.signalSensors->findSensor(GT);

  Vector<ACEplus_Vertex*> newVertices;
  if(rightBegin-ofEnd>model.maxDeNovoExonLen)
    rightBegin=ofEnd+model.maxDeNovoExonLen;
  findMateSignals(*sensor,Interval(ofEnd,rightBegin),newVertices);
  
  // Link the new vertices into the graph
  for(Vector<ACEplus_Vertex*>::iterator cur=newVertices.begin(), end=
	newVertices.end() ; cur!=end ; ++cur) {
    LightVertex *newVertex=*cur;
    //int newBegin, newEnd;
    //getContextWindow(newVertex,newBegin,newEnd);

    // Link this new vertex to the de novo vertex, to create a cryptic exon
#ifdef SANITY_CHECKS
    if(ofVertex->getType()==newVertex->getType()) INTERNAL_ERROR;
#endif
    ACEplus_Edge *edge=linkVertices(ofVertex,newVertex);
    if(edge) edge->getChange().crypticExon=true;
    
    // Link this new vertex right to the nearest annotated vertex
#ifdef SANITY_CHECKS
    if(newVertex->getType()==right->getType()) INTERNAL_ERROR;
#endif
    linkVertices(newVertex,right);
  }
}



void GraphBuilder::handleRegulatoryChanges()
{
  // First, find variants
  if(variants.isEmpty()) findVariants(variants);
  if(variants.size()==0) return;
  Vector<ExonEdge> exons;
  getAnnotatedExons(exons);

  // Variants that weaken an exon's definition, causing it to be skipped
  handleExonWeakening(variants,exons);

  // Variants that shorten an exon by weakening exon definition just inside
  // the exon
  handleExonShortening(variants,exons);

  // Variants that strengthen exon definition close to an exon, causing the
  // exon to be extended to a nearby cryptic site
  handleExonExtension(variants,exons);

  // Variants that strengthen exon definition deep within an intron, resulting
  // in a cryptic exon
  handleDeepIntronic(variants,exons);
}



double GraphBuilder::exonDefChange(Interval interval)
{
  const ContentSensor *sensor=model.contentSensors->getSensor(EXON);
  double altScore=
    model.contentSensors->score(EXON,interval.getBegin(),interval.getEnd());
  int refBegin=altToRef.mapApproximate(interval.getBegin());
  int refEnd=altToRef.mapApproximate(interval.getEnd());
  //cout<<refBegin<<" "<<refEnd<<" "<<interval<<endl;
  double refScore=sensor->scoreSubsequence(refSeq,refSeqStr,refBegin,
					   refEnd-refBegin,0);
  double altNormalized=altScore/interval.length();
  double refNormalized=refScore/(refEnd-refBegin);
  double change=altNormalized-refNormalized;
  return exp(change); // likelihood ratio (not logarithmic)
}



double GraphBuilder::exonIntronRatio(LightVertex *from,
				     LightVertex *to)
{
  int fromBegin, fromEnd, toBegin, toEnd;
  getContextWindow(from,fromBegin,fromEnd);
  getContextWindow(to,toBegin,toEnd);
  return exonIntronRatio(Interval(fromEnd,toBegin));
}



double GraphBuilder::exonIntronRatio(Interval interval)
{
  double exonScore=
    model.contentSensors->score(EXON,interval.getBegin(),interval.getEnd());
  double intronScore=
    model.contentSensors->score(INTRON,interval.getBegin(),interval.getEnd());
  double ratio=exonScore-intronScore;
  return exp(ratio); // likelihood ratio (not logarithmic)
}



void GraphBuilder::getAnnotatedExons(Vector<ExonEdge> &exons)
{
  const int numEdges=G->getNumEdges();
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *edge=G->getEdge(i);
    if(edge->isExon() && edge->isAnnotated()) {
      Interval interval(edge->getLeft()->getEnd(),
			edge->getRight()->getBegin());
      ExonEdge exonEdge;
      exonEdge.interval=interval;
      exonEdge.edge=edge;
      exons.push_back(exonEdge);
    }
  }
}



void GraphBuilder::getVariantsInExons(const Vector<Interval> &variants,
				      const Vector<ExonEdge> &exons,
				      const Vector< pair<Interval,ExonEdge> >
				      &into) 
{
  for(Vector<Interval>::const_iterator cur=variants.begin(), end=variants.end()
	; cur!=end ; ++cur) {
    const Interval &variant=*cur;
    for(Vector<ExonEdge>::const_iterator cur=exons.begin(), end=exons.end() ;
	cur!=end ; ++cur) {
      const ExonEdge &exon=*cur;
      if(exon.interval.overlaps(variant)) {
	into.push_back(pair<Interval,ExonEdge>(variant,exon));
	break;
      }
    }
  }
}



void GraphBuilder::handleExonWeakening(const Vector<Interval> &variants,
				       Vector<ExonEdge> &exons)
{
  // For each variant in an exon, evaluate change in exon definition within
  // a small window around the variant (left-5 to right+5)
  Vector< pair<Interval,ExonEdge> > variantsInExons;
  getVariantsInExons(variants,exons,variantsInExons);
  for(Vector< pair<Interval,ExonEdge> >::const_iterator 
	cur=variantsInExons.begin(), end=variantsInExons.end() ;
      cur!=end ; ++cur) {
    const pair<Interval,ExonEdge> &p=*cur;
    Interval variant=p.first;
    ExonEdge exon=p.second;
    int begin=variant.getBegin()-5;
    int end=variant.getEnd()+5;
    if(begin<exon.interval.getBegin()) begin=exon.interval.getBegin();
    if(end>exon.interval.getEnd()) end=exon.interval.getEnd();
    if(end<=begin) continue;
    //cout<<"begin="<<begin<<" end="<<end<<" exon end="<<exon.interval.getEnd()<<endl;
    const double LR=exonDefChange(Interval(begin,end));

    // If the variant substantially weakens exon definition, propose
    // skipping the exon:
    if(LR<model.EXON_WEAKENING_THRESHOLD) addExonSkippingEdge(exon);
  }
}



void GraphBuilder::getAllVertices(int from,int to,SignalType type,
				  Vector<int> &into)
{
  for(int i=from ; i<=to ; ++i) {
    LightVertex *v=G->getVertex(i);
    if(v->getType()==type) into.push_back(v->getID());
  }
}
				      


/*
  This function is only intended to be used for variants that weaken
  the exon definition score for an exon.
 */
void GraphBuilder::addExonSkippingEdge(ExonEdge exon)
{
  // Re-sort the graph and set IDs so the IDs will be stable
  G->sort();

  // First, find annotated vertices for this exon and the left and right
  // exons
  LightVertex *exonBeginVertex=exon.edge->getLeft();
  LightVertex *exonEndVertex=exon.edge->getRight();
  if(!exonBeginVertex->isSpliceSite() || 
     !exonEndVertex->isSpliceSite()) return;
  LightVertex *leftAnnotated=
    GraphBuilder::findAnnotatedVertexLeftOf(exonBeginVertex->getID());
  LightVertex *rightAnnotated=
    GraphBuilder::findAnnotatedVertexRightOf(exonEndVertex->getID());

  // Get all donor sites back to the annotated site on the left, and all
  // acceptor sites up to the annotated site on the right
  Vector<int> donors, acceptors;
  getAllVertices(leftAnnotated->getID(),exonBeginVertex->getID(),GT,donors);
  getAllVertices(exonEndVertex->getID(),rightAnnotated->getID(),AG,acceptors);
  
  // Link all pairs of donors/acceptors immediately enclosing this exon
  for(Vector<int>::iterator cur=donors.begin(), end=donors.end() ;
      cur!=end ; ++cur) {
    LightVertex *donor=G->getVertex(*cur);
    for(Vector<int>::iterator cur=acceptors.begin(), end=acceptors.end() ;
	cur!=end ; ++cur) {
      LightVertex *acceptor=G->getVertex(*cur);
#ifdef SANITY_CHECKS
      if(donor->getType()==acceptor->getType()) INTERNAL_ERROR;
#endif
      ACEplus_Edge *edge=linkVertices(donor,acceptor);
      if(edge) {
	edge->getChange().exonSkipping=true;
	edge->getChange().regulatoryChange=true; }
    }    
  }
}



bool GraphBuilder::isInExon(Interval variant,Vector<ExonEdge> &exons)
{
  for(Vector<ExonEdge>::iterator cur=exons.begin(), end=exons.end() ;
      cur!=end ; ++cur) {
    const ExonEdge &exon=*cur;
    if(exon.interval.overlaps(variant)) return true;
  }
  return false;
}



void GraphBuilder::handleExonExtension(const Vector<Interval> &variants,
				       Vector<ExonEdge> &exons)
{
  // Consider each variant close to an annotated exon
  for(Vector<Interval>::const_iterator cur=variants.begin(), end=variants.end()
	; cur!=end ; ++cur) {
    const Interval &variant=*cur;
    if(isInExon(variant,exons)) continue;
    for(Vector<ExonEdge>::iterator cur=exons.begin(), end=exons.end() ;
	cur!=end ; ++cur) {
      const ExonEdge &exon=*cur;
      if(variant.distanceTo(exon.interval)<=model.MAX_SPLICE_SHIFT) {
	int begin=variant.getBegin()-5;
	int end=variant.getEnd()+5;
	if(begin<exon.interval.getEnd()) {
	  if(end>exon.interval.getBegin()) end=exon.interval.getBegin(); }
	else
	  if(begin<exon.interval.getEnd()) begin=exon.interval.getEnd();
        if(end<=begin) continue;
	const double LR=exonDefChange(Interval(begin,end));
	
	// If the variant substantially strengthens exon definition beyond
	// the end of the exon, propose lengthening the exon to a cryptic site
	if(LR>model.EXON_STRENGTHENING_THRESHOLD) lengthenExon(variant,exon);
      }
    }
  }
}



void GraphBuilder::lengthenExon(const Interval &variant,const ExonEdge &exon)
{
  if(variant.getEnd()<exon.interval.getBegin())
    lengthenExonLeft(variant,exon);
  else if(variant.getBegin()>exon.interval.getEnd())
    lengthenExonRight(variant,exon);
  else INTERNAL_ERROR;
}



void GraphBuilder::scan(const Interval &scanWindow,SignalType signalType,
			Vector<ACEplus_Vertex*> &into)
{
  SignalSensor &sensor=*model.signalSensors->findSensor(signalType);
  const int consensusOffset=sensor.getConsensusOffset();
  const int consensusLen=sensor.getConsensusLength();
  const String &substrate=projected.getSubstrate();
  const Strand strand=projected.getStrand();
  const int sensorLen=sensor.getContextWindowLength();
  for(int pos=scanWindow.getBegin() ; pos<scanWindow.getEnd() ; ++pos) {
    if(pos+sensorLen>altSeq.getLength()) continue;
    if(sensor.consensusOccursAt(altSeqStr,pos+consensusOffset)) {
      const double altScore=sensor.getLogP(altSeq,altSeqStr,pos);
      if(altScore<sensor.getCutoff()) continue;
      ACEplus_Vertex *v=newVertex(substrate,signalType,pos+consensusOffset,
			       pos+consensusOffset+consensusLen,
			       altScore,strand,G->getNumVertices());
      if(!v) continue;
      G->addVertex(v);
      into.push_back(v);
      v->setThreshold(sensor.getCutoff());
    }
  }
}



void GraphBuilder::lengthenExonLeft(const Interval &variant,
				    const ExonEdge &exon)
{
  // Scan for new exon ends
  LightVertex *exonBegin=exon.edge->getLeft();
  int scanBegin=exon.interval.getBegin()-model.MAX_SPLICE_SHIFT;
  int scanEnd=variant.getBegin();
  if(scanBegin>=scanEnd) return;
  Vector<ACEplus_Vertex*> found;
  scan(Interval(scanBegin,scanEnd),AG,found);

  // Link the new ends into the graph
  for(Vector<ACEplus_Vertex*>::iterator cur=found.begin(), end=found.end() ;
      cur!=end ; ++cur) {
    LightVertex *newVertex=*cur;

    // Add exon edges
    for(Vector<LightEdge*>::iterator cur=exonBegin->getEdgesOut().begin(),
	  end=exonBegin->getEdgesOut().end() ; cur!=end ; ++cur) {
      LightVertex *linkTo=(*cur)->getRight();
#ifdef SANITY_CHECKS
      if(newVertex->getType()==linkTo->getType()) INTERNAL_ERROR;
#endif
      ACEplus_Edge *edge=linkVertices(newVertex,linkTo);
      if(edge) edge->getChange().regulatoryChange=true;}

    // Add intron edges
    for(Vector<LightEdge*>::iterator cur=exonBegin->getEdgesIn().begin(),
	  end=exonBegin->getEdgesIn().end() ; cur!=end ; ++cur) {
      LightVertex *linkTo=(*cur)->getLeft();
#ifdef SANITY_CHECKS
      if(linkTo->getType()==newVertex->getType()) INTERNAL_ERROR;
#endif
      ACEplus_Edge *edge=linkVertices(linkTo,newVertex);
      if(edge) edge->getChange().regulatoryChange=true;}
  }
}



void GraphBuilder::lengthenExonRight(const Interval &variant,
				     const ExonEdge &exon)
{
  // Scan for new exon ends
  LightVertex *exonEnd=exon.edge->getRight();
  int scanBegin=variant.getEnd();
  int scanEnd=exon.interval.getEnd()+model.MAX_SPLICE_SHIFT;
  if(scanBegin>=scanEnd) return;
  Vector<ACEplus_Vertex*> found;
  scan(Interval(scanBegin,scanEnd),GT,found);

  // Link the new ends into the graph
  for(Vector<ACEplus_Vertex*>::iterator cur=found.begin(), end=found.end() ;
      cur!=end ; ++cur) {
    LightVertex *newVertex=*cur;

    // Add intron edges
    for(Vector<LightEdge*>::iterator cur=exonEnd->getEdgesOut().begin(),
	  end=exonEnd->getEdgesOut().end() ; cur!=end ; ++cur) {
      LightVertex *linkTo=(*cur)->getRight();
#ifdef SANITY_CHECKS
      if(newVertex->getType()==linkTo->getType()) INTERNAL_ERROR;
#endif
      ACEplus_Edge *edge=linkVertices(newVertex,linkTo);
      if(edge) edge->getChange().regulatoryChange=true;}

    // Add exon edges
    for(Vector<LightEdge*>::iterator cur=exonEnd->getEdgesIn().begin(),
	  end=exonEnd->getEdgesIn().end() ; cur!=end ; ++cur) {
      LightVertex *linkTo=(*cur)->getLeft();
#ifdef SANITY_CHECKS
      if(linkTo->getType()==newVertex->getType()) INTERNAL_ERROR;
#endif
      ACEplus_Edge *edge=linkVertices(linkTo,newVertex);
      if(edge) edge->getChange().regulatoryChange=true;}
  }
}



int GraphBuilder::distanceOfExonicVariantToEnd(const Interval &variant,
					       const ExonEdge &exon)
{
  const int dist1=exon.interval.getEnd()-variant.getBegin();
  const int dist2=variant.getEnd()-exon.interval.getBegin();
  const int minDist=dist1<dist2 ? dist1 : dist2;
#ifdef SANITY_CHECKS
  if(minDist<0) INTERNAL_ERROR;
#endif
  return minDist;
}



void GraphBuilder::handleExonShortening(const Vector<Interval> &variants,
					Vector<ExonEdge> &exons)
{
  // Consider each variant within an annotated exon, close to the end
  for(Vector<Interval>::const_iterator cur=variants.begin(), end=variants.end()
	; cur!=end ; ++cur) {
    const Interval &variant=*cur;
    for(Vector<ExonEdge>::iterator cur=exons.begin(), end=exons.end() ;
	cur!=end ; ++cur) {
      const ExonEdge &exon=*cur;
      if(!exon.interval.overlaps(variant)) continue;
      if(distanceOfExonicVariantToEnd(variant,exon)<=model.MAX_SPLICE_SHIFT) {
	int begin=variant.getBegin()-5;
	int end=variant.getEnd()+5;
	if(begin<exon.interval.getBegin()) begin=exon.interval.getBegin();
	if(end>exon.interval.getEnd()) end=exon.interval.getEnd();
        if(end<=begin) continue;
	const double LR=exonDefChange(Interval(begin,end));

	// If the variant substantially weakens exon definition, propose
	// shortening the exon to a cryptic splice site within the exon
	if(LR<model.EXON_WEAKENING_THRESHOLD) shortenExon(variant,exon);
      }
    }
  }
}



void GraphBuilder::shortenExon(const Interval &variant,const ExonEdge &exon)
{
  const int distToRightEnd=exon.interval.getEnd()-variant.getBegin();
  const int distToLeftEnd=variant.getEnd()-exon.interval.getBegin();
  if(distToRightEnd<=distToLeftEnd)
    shortenExonFromTheRight(variant,exon);
  else 
    shortenExonFromTheLeft(variant,exon);
}



void GraphBuilder::shortenExonFromTheRight(const Interval &variant,
					   const ExonEdge &exon)
{
  // Scan for new exon ends
  LightVertex *exonEnd=exon.edge->getRight();
  int scanBegin=exon.interval.getEnd()-model.MAX_SPLICE_SHIFT;
  int scanEnd=variant.getBegin();
  if(scanBegin>=scanEnd) return;
  Vector<ACEplus_Vertex*> found;
  scan(Interval(scanBegin,scanEnd),GT,found);

  // Link the new ends into the graph
  for(Vector<ACEplus_Vertex*>::iterator cur=found.begin(), end=found.end() ;
      cur!=end ; ++cur) {
    LightVertex *newVertex=*cur;

    // Add intron edges
    for(Vector<LightEdge*>::iterator cur=exonEnd->getEdgesOut().begin(),
	  end=exonEnd->getEdgesOut().end() ; cur!=end ; ++cur) {
      LightVertex *linkTo=(*cur)->getRight();
#ifdef SANITY_CHECKS
      if(newVertex->getType()==linkTo->getType()) INTERNAL_ERROR;
#endif
      if(newVertex->getEnd()>=linkTo->getBegin()) continue;
      ACEplus_Edge *edge=linkVertices(newVertex,linkTo);
      if(edge) edge->getChange().regulatoryChange=true;}

    // Add exon edges
    for(Vector<LightEdge*>::iterator cur=exonEnd->getEdgesIn().begin(),
	  end=exonEnd->getEdgesIn().end() ; cur!=end ; ++cur) {
      LightVertex *linkTo=(*cur)->getLeft();
#ifdef SANITY_CHECKS
      if(linkTo->getType()==newVertex->getType()) INTERNAL_ERROR;
#endif
      if(linkTo->getEnd()>=newVertex->getBegin()) continue;
      ACEplus_Edge *edge=linkVertices(linkTo,newVertex);
      if(edge) edge->getChange().regulatoryChange=true;}
  }
}



void GraphBuilder::shortenExonFromTheLeft(const Interval &variant,
					  const ExonEdge &exon)
{
  // Scan for new exon ends
  LightVertex *exonBegin=exon.edge->getLeft();
  int scanBegin=variant.getEnd();
  int scanEnd=exon.interval.getBegin()+model.MAX_SPLICE_SHIFT;
  if(scanBegin>=scanEnd) return;
  Vector<ACEplus_Vertex*> found;
  scan(Interval(scanBegin,scanEnd),AG,found);

  // Link the new ends into the graph
  for(Vector<ACEplus_Vertex*>::iterator cur=found.begin(), end=found.end() ;
      cur!=end ; ++cur) {
    LightVertex *newVertex=*cur;

    // Add exon edges
    for(Vector<LightEdge*>::iterator cur=exonBegin->getEdgesOut().begin(),
	  end=exonBegin->getEdgesOut().end() ; cur!=end ; ++cur) {
      LightVertex *linkTo=(*cur)->getRight();
#ifdef SANITY_CHECKS
      if(newVertex->getType()==linkTo->getType()) INTERNAL_ERROR;
#endif
      if(newVertex->getEnd()>=linkTo->getBegin()) continue;
      ACEplus_Edge *edge=linkVertices(newVertex,linkTo);
      if(edge) edge->getChange().regulatoryChange=true;}

    // Add intron edges
    for(Vector<LightEdge*>::iterator cur=exonBegin->getEdgesIn().begin(),
	  end=exonBegin->getEdgesIn().end() ; cur!=end ; ++cur) {
      LightVertex *linkTo=(*cur)->getLeft();
#ifdef SANITY_CHECKS
      if(linkTo->getType()==newVertex->getType()) INTERNAL_ERROR;
#endif
      if(linkTo->getEnd()>=newVertex->getBegin()) continue;
      ACEplus_Edge *edge=linkVertices(linkTo,newVertex);
      if(edge) edge->getChange().regulatoryChange=true;}
  }
}



void GraphBuilder::handleDeepIntronic(const Vector<Interval> &variants,
				      Vector<ExonEdge> &exons)
{
  for(Vector<Interval>::const_iterator cur=variants.begin(), end=variants.end()
	; cur!=end ; ++cur) {
    const Interval &variant=*cur;
    if(isInExon(variant,exons)) continue;
    const int numExons=exons.size();
    for(int i=0 ; i<numExons-1 ; ++i) {
      const ExonEdge &thisExon=exons[i], &nextExon=exons[i+1];
      Interval intron(thisExon.interval.getEnd(),nextExon.interval.getBegin());
      if(intron.overlaps(variant)) {
	const int dist1=variant.getBegin()-intron.getBegin();
	const int dist2=intron.getEnd()-variant.getEnd();
	if(dist1<model.MIN_INTRON_LEN || dist2<dist1<model.MIN_INTRON_LEN)
	  continue;
	int begin=variant.getBegin()-5;
	int end=variant.getEnd()+5;
	if(end<begin || begin<intron.getBegin() || end>intron.getEnd())
	  continue;
	const double LR=exonDefChange(Interval(begin,end));

	// The variant must result in a substantial increase in exon
	// definition score, relative to the reference, for us to predict
	// a cryptic exon:
	if(LR>model.EXON_STRENGTHENING_THRESHOLD)
	  proposeCrypticExons(variant,thisExon,nextExon);
      }
    }
  }
}



void GraphBuilder::proposeCrypticExons(const Interval &variant,
				       const ExonEdge &prevExon,
				       const ExonEdge &nextExon)
{
  // Scan for splice sites on either side of the variant
  int scanBegin=variant.getBegin()-model.maxDeNovoExonLen;
  int scanEnd=variant.getEnd()+model.maxDeNovoExonLen;
  Vector<ACEplus_Vertex*> acceptors, donors;
  scan(Interval(scanBegin,variant.getBegin()),AG,acceptors);
  scan(Interval(variant.getEnd(),scanEnd),GT,donors);
  
  // Link splice sites into the graph
  int prevExonEndID=prevExon.edge->getRight()->getID();
  int nextExonBeginID=nextExon.edge->getLeft()->getID();
  for(Vector<ACEplus_Vertex*>::iterator cur=acceptors.begin(), 
	end=acceptors.end() ; cur!=end ; ++cur) {
    LightVertex *acceptor=*cur;
    for(Vector<ACEplus_Vertex*>::iterator cur=donors.begin(), 
	  end=donors.end() ; cur!=end ; ++cur) {
      LightVertex *donor=*cur;
      if(donor->getBegin()-acceptor->getEnd()>model.maxDeNovoExonLen)
	continue;

      // Apply a filter on the exon definition score for this interval
      if(exonIntronRatio(acceptor,donor)<model.MIN_EXON_INTRON_RATIO)
	continue;

      // Create the edge denoting the cryptic exon
#ifdef SANITY_CHECKS
      if(acceptor->getType()==donor->getType()) INTERNAL_ERROR;
#endif
      ACEplus_Edge *edge=linkVertices(acceptor,donor);
      if(edge) {
	edge->getChange().crypticExon=true;
	edge->getChange().regulatoryChange=true;}

      // Create edges for the introns left and right of the cryptic exon
      Vector<int> leftTargets, rightTargets;
      getAllVertices(prevExonEndID,nextExonBeginID,GT,leftTargets);
      getAllVertices(prevExonEndID,nextExonBeginID,AG,rightTargets);
      for(Vector<int>::const_iterator cur=leftTargets.begin(), end=
	    leftTargets.end() ; cur!=end ; ++cur) {
	LightVertex *v=G->getVertex(*cur);
#ifdef SANITY_CHECKS
	if(v->getType()==acceptor->getType()) INTERNAL_ERROR;
#endif
	if(acceptor->getBegin()-v->getEnd()>=model.MIN_INTRON_LEN)
	  linkVertices(v,acceptor); }
      for(Vector<int>::const_iterator cur=rightTargets.begin(), end=
	    rightTargets.end() ; cur!=end ; ++cur) {
	LightVertex *v=G->getVertex(*cur);
#ifdef SANITY_CHECKS
	if(donor->getType()==v->getType()) INTERNAL_ERROR;
#endif
	if(v->getBegin()-donor->getEnd()>=model.MIN_INTRON_LEN)
	  linkVertices(donor,v); }
    }
  }
}



bool GraphBuilder::buildGraph(bool strict)
{
  // First, build a basic graph from the projected annotation
  if(!buildTranscriptGraph()) return false;
  //cout<<"GRAPH 1\n"<<*G<<endl;
  cout<<G->getNumVertices()<<" vertices in graph"<<endl;

  if(!strict) {
    // Add vertices & edges to address broken structures
    cout<<"handling broken sites"<<endl;
    if(signals.anyBroken()) handleBrokenSites();
    cout<<G->getNumVertices()<<" vertices in graph"<<endl;

    // Add novel sites created by genetic variants
    cout<<"predicting de novo sites"<<endl;
    if(model.allowDeNovoSites) handleDeNovoSites();
    cout<<G->getNumVertices()<<" vertices in graph"<<endl;

    // Add vertices/edges to accommodate structure changes due to ESE variants
    cout<<"predicting regulatory changes"<<endl;
    if(model.allowRegulatoryChanges) handleRegulatoryChanges();
    cout<<G->getNumVertices()<<" vertices in graph"<<endl;

    // Prune away any vertex or edge not reachable from both ends
    cout<<G->getNumVertices()<<" vertices before pruning"<<endl;
    //cout<<*G<<endl;
    pruneUnreachable(*G);
    cout<<G->getNumVertices()<<" vertices after pruning"<<endl;
  }

  // Mark any edges that result in intron retention
  markIntronRetentions();

  // Score edges
  cout<<"scoring edges"<<endl;
  const int numEdges=G->getNumEdges();
  cout<<numEdges<<" edges"<<endl;
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *edge=G->getEdge(i);
    scoreEdge(dynamic_cast<ACEplus_Edge*>(edge));
    if(!isFinite(edge->getScore())) {
      cout<<"-inf edge score: "<<*edge<<endl;
      G->dropEdge(i);
    }
  }
  cout<<"deleting null edges"<<endl;
  G->deleteNullEdges();

  // Determine whether mapping was exact
  cout<<"checking for any changes"<<endl;
  if(!allVerticesAreAnnotated()) changes=true;
  cout<<"done checking"<<endl;
  return true;
}



void GraphBuilder::getAnnotatedEdges(Vector<LightEdge*> &into)
{
  const int numEdges=G->getNumEdges();
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *edge=G->getEdge(i);
    if(edge->isAnnotated()) into.push_back(edge);
  }
}



bool GraphBuilder::coversAnnotatedIntron(const LightEdge &edge,
					 Vector<LightEdge*> &annotatedEdges)
{
  const Interval &I=edge.asInterval();
  for(Vector<LightEdge*>::iterator cur=annotatedEdges.begin(),
	end=annotatedEdges.end() ; cur!=end ; ++cur) {
    LightEdge *other=*cur;
    if(other->isIntron() && I.contains(other->asInterval())) return true;
    if(other->isIntron() && I.contains(other->asInterval())) {
      cout<<"XXX "<<edge<<" covers "<<*other<<endl; // ###
    }
  }
  return false;
}


void GraphBuilder::markIntronRetentions()
{
  Vector<LightEdge*> annotated;
  getAnnotatedEdges(annotated);
  const int numEdges=G->getNumEdges();
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *edge=G->getEdge(i);
    if(edge->isExon() && coversAnnotatedIntron(*edge,annotated))
      dynamic_cast<ACEplus_Edge*>(edge)->getChange().intronRetention=true;
  }
}




