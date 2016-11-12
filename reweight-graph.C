/****************************************************************
 reweight-graph.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "BOOM/Constants.H"
#include "BOOM/CommandLine.H"
#include "BOOM/WigBinary.H"
#include "BOOM/Array1D.H"
#include "BOOM/PriorityTree.H"
#include "RnaJunctions.H"
#include "LightGraph.H"
#include "IntronDepthProfile.H"
using namespace std;
using namespace BOOM;


const float PSEUDOCOUNT=1;

class Application {
public:
  void main(int argc,char *argv[]);
protected:
  float dnaWeight, rnaWeight, scalingFactor, maxDepth;
  int L; // input sequence length
  int V; // number of vertices
  int E; // number of edges
  bool exonSpliceFlag; // exons are marked as supported only if their
                       // splice sites are supported (no depth requirement)
  ScaledWigBinary *wig;
  RnaJunctions *junctions;
  float maxJunctionDepth, wigMax, wigAve;
  LightGraph *G;
  IntronDepthProfile *intronProfile;
  Array1D<float> PSA, negPSA;
  float computePSA();
  void processGraph();
  void processEdges();
  void processVertices();
  void removeNontermIntergenic();
  void markSupportedExons();
  void markSupportedExons_splice();
};



int main(int argc,char *argv[]) {
  try {	Application app;  app.main(argc,argv); }
  catch(const char *p) { cerr << p << endl; return -1; }
  catch(const string &msg) { cerr << msg.c_str() << endl; return -1; }
  catch(const exception &e) {
    cerr << "STL exception caught in main:\n" << e.what() << endl;
    return -1;
  }
  catch(...) {
    cerr << "Unknown exception caught in main" << endl;
    return -1;
  }
  return 0;
}



void Application::main(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"ie");
  if(cmd.numArgs()!=6)
    throw BOOM::String("\n\
reweight-graph [-i] <in.graph> <*.pileup> <*.junctions> <genezilla-weight> <rna-weight> <out.graph>\n\
   -i = remove nonterminal intergenic edges\n\
   -e = base exon support flag on splicing support\n\
");
  String inGraphFile=cmd.arg(0);
  String pileupFile=cmd.arg(1);
  String junctionFile=cmd.arg(2);
  dnaWeight=cmd.arg(3).asFloat();
  rnaWeight=cmd.arg(4).asFloat();
  String outfile=cmd.arg(5);
  bool dashI=cmd.option('i');
  exonSpliceFlag=cmd.option('e');

  // Load data
  G=new LightGraph(inGraphFile);
  if(dashI) removeNontermIntergenic();
  wig=new ScaledWigBinary(pileupFile);
  wig->useLog();
  L=wig->getLength();
  junctions=new RnaJunctions(junctionFile);
  junctions->logify(); // ###
  intronProfile=new IntronDepthProfile(*junctions,L);

  // Compute some global things
  maxJunctionDepth=intronProfile->getMax();//junctions->getMaxDepth();
  float wigMin;
  wig->getExtrema(wigMin,wigMax);
  maxDepth=max(wigMax,maxJunctionDepth);
  scalingFactor=1/(wigMax+PSEUDOCOUNT);
  wigAve=wig->ave();
  wig->changeScale(scalingFactor);
  wigAve*=scalingFactor;
  intronProfile->rescale(1/(maxJunctionDepth+1));//###
  V=G->getNumVertices();
  E=G->getNumEdges();
  computePSA();

  // Process the graph
  processGraph();

  G->save(outfile);
}



void Application::processGraph()
{
  processEdges();
  processVertices();
  if(exonSpliceFlag) markSupportedExons_splice();
  else markSupportedExons();
}



void Application::processEdges()
{
  // Process all elements in the graph
  for(int i=0 ; i<E ; ++i) {
    LightEdge *edge=G->getEdge(i);
    if(!edge) continue;
    int begin=edge->getBegin(), end=edge->getEnd();
    if(begin<0) begin=0;
    if(end>L) end=L;
    int len=edge->getLength();
    if(edge->isCoding()) {
      edge->subsumeVertexScores();
      float rnaScore=PSA[end-1]-(begin>0 ? PSA[begin-1] : 0);
      Array1D<int> frames=edge->getFrames();
      for(int i=0 ; i<frames.size() ; ++i) {
	float dnaScore=edge->getScore(frames[i]);
	float newScore=dnaWeight*dnaScore+rnaWeight*rnaScore;
	edge->setScore(frames[i],newScore);
      }
    }
    else if(edge->isIntron()) {
      float depth=junctions->getDepth(begin,end);
      edge->setSupport(depth>0);
      float rnaScore=log((depth+PSEUDOCOUNT)/maxJunctionDepth)*len;
      Array1D<int> frames=edge->getFrames();
      for(int i=0 ; i<frames.size() ; ++i) {
	float dnaScore=edge->getScore(frames[i]);
	float newScore=dnaWeight*dnaScore+rnaWeight*rnaScore;
	edge->setScore(frames[i],newScore);
      }      
    }
    else { // INTERGENIC
      float rnaScore=negPSA[end-1]-(begin>0 ? negPSA[begin-1] : 0);
      float dnaScore=edge->getScore(0);
      float newScore=dnaWeight*dnaScore+rnaWeight*rnaScore;
      edge->setScore(0,newScore);
    }
  }
}



void Application::processVertices()
{
  for(int i=0 ; i<V ; ++i) {
    LightVertex *vertex=G->getVertex(i);
    if(!vertex) continue;
    SignalType t=vertex->getType();
    if(beginsCoding(t) || endsCoding(t))
      vertex->setScore(0);
  }
}



float Application::computePSA()
{
  PSA.resize(L);
  negPSA.resize(L);
  static const float pseudocount=1e-8;
  PSA[0]=log(wig->read(0)+pseudocount);
  negPSA[0]=log((1-wig->read(0))*(1-(*intronProfile)[0])+pseudocount);
  for(int i=1 ; i<L ; ++i) {
    float wigVal=wig->read(i);
    PSA[i]=PSA[i-1]+log(wigVal+pseudocount);
    negPSA[i]=
      negPSA[i-1]+log((1-wigVal)*(1-(*intronProfile)[i])+pseudocount);
  }
}



void Application::removeNontermIntergenic()
{
  int N=G->getNumEdges();
  for(int i=0 ; i<N ; ++i) {
    LightEdge *edge=G->getEdge(i);
    if(!edge->isIntergenic()) continue;
    LightVertex *left=edge->getLeft();
    LightVertex *right=edge->getRight();
    if(left->getEdgesIn().isEmpty() ||
       right->getEdgesOut().isEmpty()) continue;
    Array1D<int> frames=edge->getFrames();
    int n=frames.size();
    for(int j=0 ; j<n ; ++j) edge->setScore(frames[j],NEGATIVE_INFINITY);
  }
}



void Application::markSupportedExons()
{
  const int L=wig->getLength();
  int numV=G->getNumVertices(), nextVertexInd=0;
  if(numV==0) return;
  int nextVertexIndex=0;
  EdgeEndComparator edgeEndCmp;
  PriorityTree<LightEdge*> heap(edgeEndCmp);
  heap.enableDuplicates();
  for(int currentPos=0 ; currentPos<L ; ++currentPos) {
    float reads=wig->read(currentPos);
    if(reads==0) {
      heap.purge();
    }
    else {
      // Pass vertices as necessary
      while(G->getVertex(nextVertexIndex)->getEnd()+1<currentPos) 
	++nextVertexIndex;

      // Add new edges to the heap
      for(int i=nextVertexIndex ; i<numV ; ++i) {
	LightVertex *v=G->getVertex(i);
	if(v->getBegin()>currentPos) break;
	if(v->getEnd()+1==currentPos) {
	  Vector<LightEdge*> &edgesOut=v->getEdgesOut();
	  for(Vector<LightEdge*>::iterator cur=edgesOut.begin(),
		end=edgesOut.end() ; cur!=end ; ++cur) {
	    heap.insert(*cur);
	  }
	}
      }

      // Graduate fully covered edges from the heap
      while(!heap.isEmpty()) {
	LightEdge *edge=heap.peekMin();
	if(edge->getEnd()<=currentPos+1) {
	  heap.extractMin();
	  edge->setSupport(true);
	}
	else break;
      }
    }
  }
}



void Application::markSupportedExons_splice()
{
  markSupportedExons();
  int numEdges=G->getNumEdges();
  for(int i=0 ; i<numEdges ; ++i) {
    LightEdge *edge=G->getEdge(i);
    if(!edge->isSupported()) continue;
    if(!edge->isCoding()) continue;
    LightVertex *left=edge->getLeft(), *right=edge->getRight();
    SignalType leftType=dropStrand(left->getType());
    SignalType rightType=dropStrand(right->getType());
    bool supported=true;
    if(leftType==GT || leftType==AG && 
       junctions->getSpliceInDepth(left->getBegin()+2)==0) supported=false;
    if(rightType==GT || rightType==AG && 
       junctions->getSpliceInDepth(right->getBegin())==0) supported=false;
    //if(dropStrand(edge->getType())==SINGLE_EXON) supported=true;
    edge->setSupport(supported);
  }
}
