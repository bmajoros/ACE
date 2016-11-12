/****************************************************************
 Get-Classifier-Examples
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "genezilla.H"
#include "GeneZilla.H"
#include "EdgeFactory.H"
#include "GarbageCollector.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/Random.H"
#include "BOOM/Constants.H"
#include "NthOrderStringIterator.H"
#include <assert.h>
using namespace std;

/****************************************************************
                         defines and pragmas
 */
#ifndef EXPLICIT_GRAPHS
#error Please edit the file genezilla.H, uncomment the line defining EXPLICIT_GRAPHS, issue a "make clean" and then recompile this project
#endif
const int LARGE_PRIME=3121;

/****************************************************************
                               globals
 */
static const char *PROGRAM_NAME="get-classifier-examples";
static const char *VERSION="1.0";
Alphabet alphabet;
int frame; // ### CAUTION: this is required by older code; to be removed
GarbageCollector garbageCollector;


/****************************************************************
                             IntervalSet
 */
typedef pair<int,int> Interval;
int hashInterval(const Interval &i) {return i.first*i.second;}
struct IntervalSetCmp
{
  bool operator()(const Interval &a,const Interval &b) const
  {
    //return hashInterval(a)<hashInterval(b);
    return a.first<b.first || a.first==b.first && a.second<b.second;
  }
};
typedef BOOM::Set<Interval,IntervalSetCmp> IntervalSet;


/****************************************************************
                           class Application
 */
class Application
{
  BOOM::StringMap<double> hexamerLogRatios;
  ofstream osInitial, osInternal, osFinal, osSingle;

  void output(ContentType,double lengthProb,double firstScore,
	      double secondScore,double hmmScore,int category);
  int countExons(TranscriptList &);
  void getIntervals(TranscriptList &,IntervalSet &);
  void outputPositives(TranscriptList &,ostream &,GeneZilla &,
		       const BOOM::String &seq,ParseGraph &,
		       const Sequence &);
  void outputPositive(BOOM::GffExon &,ostream &,GeneZilla &,
		      const BOOM::String &seq,ParseGraph &,
		      const Sequence &);
  ContentType exonTypeToContentType(ExonType,Strand);
  void exonCoordsToSigPos(ContentType,int begin,int end,int &pos1,
			  int &pos2);
  void sigPosToExonCoords(ContentType,int &begin,int &end,int pos1,
			  int pos2);
  double scoreSignal(SignalType signalType,int pos,ParseGraph &graph,
		     const Sequence &,const BOOM::String &seqStr,
		     GeneZilla &);
  int getBackgrounds(const BOOM::String &filename,
		     BOOM::StringMap<int> &backgroundCounts);
  void countHexamers(int begin,int end,const BOOM::String &seqStr,
		     Strand,BOOM::StringMap<int> &counts);
  void countHexamers(const char *,int len,BOOM::StringMap<int> &counts);
  double scoreHexamers(int begin,int end,const BOOM::String &seqStr,
		       Strand);
  double scoreHMM(int begin,int end,const BOOM::String &seqStr,
		  const Sequence &,
		  ContentType,GeneZilla &);
  void getLogRatios(const BOOM::String &backgroundFilename,
		    BOOM::Map<BOOM::String,TranscriptList*> &gff,
		    const BOOM::String &substrateFilename);
  void outputNegatives(TranscriptList &knownGenes,ostream &os,
		       ParseGraph &orfGraph,const BOOM::String &seqString,
		       const Sequence &seq,GeneZilla &);
public:
  Application();
  void main(int argc,char *argv[]);
};



/****************************************************************
                             main()
 */
int main(int argc,char *argv[])
  {
    try
      {
	randomize();
	Application app;
	app.main(argc,argv);
	return 0;
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application() 
  : hexamerLogRatios(LARGE_PRIME),
    osInitial("initial-exons.1"),
    osInternal("internal-exons.1"),
    osFinal("final-exons.1"),
    osSingle("single-exons.1")
{
}



void Application::main(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw BOOM::String(PROGRAM_NAME)+" <exons.gff> <substrates.fasta> <background-dna.fasta> <*.iso>";
  BOOM::String gffFilename=cmd.arg(0);
  BOOM::String contigsFilename=cmd.arg(1);
  BOOM::String backgroundFilename=cmd.arg(2);
  BOOM::String isoFilename=cmd.arg(3);

  // Load the correct gene model from the GFF file
  cerr<<"loading gene model..."<<endl;
  BOOM::GffReader gffReader(gffFilename);
  BOOM::Map<BOOM::String,TranscriptList*> &genesByContig=
    *gffReader.loadByContig();

  // Compute hexamer log likelihood ratios
  cerr<<"compute hexamer log likelihood ratios..."<<endl;
  getLogRatios(backgroundFilename,genesByContig,contigsFilename);

  // Process all the contigs
  cerr<<"processing contigs..."<<endl;
  BOOM::FastaReader fastaReader(contigsFilename);
  BOOM::String defline, seqString;
  int transcriptId=1000;
  float gcContent;
  while(fastaReader.nextSequence(defline,seqString))
    {
      Sequence seq(seqString,DnaAlphabet::global);

      // Parse the substrate ID out of the defline
      BOOM::String contigId, remainder;
      BOOM::FastaReader::parseDefline(defline,contigId,remainder);
      cerr<<"processing contig "<<contigId<<"...("
	  <<seqString.length()<<" bp)..."<<endl;

      // Get the ORF graph from GeneZilla
      cerr<<"identifying open reading frames..."<<endl;
      EdgeFactory edgeFactory;
      GeneZilla genezilla(PROGRAM_NAME,VERSION,edgeFactory,transcriptId);
      ParseGraph &orfGraph=
	genezilla.parse(isoFilename,seq,seqString,gcContent);
      orfGraph.setVertexIndices();

      // Get all the known genes on this contig from the GFF
      TranscriptList *knownGenes=genesByContig[contigId];
      if(!knownGenes) continue;

      // Output positive exon examples
continue;//###
      cerr<<"generating positive examples..."<<endl;
      outputPositives(*knownGenes,cout,genezilla,seqString,orfGraph,seq);

      // Find negative exon examples in the noncoding regions
      cerr<<"generating negative examples..."<<endl;
      outputNegatives(*knownGenes,cout,orfGraph,seqString,seq,genezilla);
    }
  
  // Clean up
  BOOM::Map<BOOM::String,TranscriptList*>:: iterator 
    cur=genesByContig.begin(), end=genesByContig.end();
  for(; cur!=end ; ++cur)
    {
      if(!(*cur).second) continue;
      TranscriptList &list=*(*cur).second;
      BOOM::Vector<BOOM::GffTranscript*>::iterator cur=list.begin(),
	end=list.end();
      for(; cur!=end ; ++cur) delete *cur;
      delete &list;
    }
  delete &genesByContig;
}



void Application::outputPositives(TranscriptList &transcripts,
				  ostream &os,GeneZilla &genezilla,
				  const BOOM::String &seqStr,
				  ParseGraph &graph,
				  const Sequence &seq)
{
  // Iterate through the transcripts
  TranscriptList::iterator tCur=transcripts.begin(), 
    tEnd=transcripts.end();
  for(; tCur!=tEnd ; ++tCur)
    {
      BOOM::GffTranscript &transcript=**tCur;

      // Iterate through the transcript's exons:
      BOOM::Vector<BOOM::GffExon*>::iterator eCur=transcript.getExons(),
	eEnd=transcript.getExonsEnd();
      for(; eCur!=eEnd ; ++eCur)
	{
	  BOOM::GffExon &exon=**eCur;
	  
	  // Output the features for this exon:
	  outputPositive(exon,os,genezilla,seqStr,graph,seq);
	}
    }
}



void Application::outputPositive(BOOM::GffExon &exon,ostream &os,
				 GeneZilla &genezilla,
				 const BOOM::String &seqStr,
				 ParseGraph &graph,const Sequence &seq)
{
  ExonType exonType=exon.getExonType();
  Strand strand=exon.getStrand();
  int begin=exon.getBegin();
  int end=exon.getEnd();
  int length=end-begin;
  if(begin<50 || end+50>seqStr.length()) return;

  // GET LENGTH PROBABILITY
  ContentType contentType=GeneZilla::exonTypeToContentType(exonType,strand);
  DiscreteDistribution &distribution=genezilla.getDistribution(contentType);
  double lengthProb=distribution.getLogP(length);
  
  // GET SIGNAL SCORES
  SignalType leftSignal=::leftSignal(contentType);
  SignalType rightSignal=::rightSignal(contentType);
  int leftPos,rightPos;
  exonCoordsToSigPos(contentType,begin,end,leftPos,rightPos);
  double firstScore=
    scoreSignal(leftSignal,leftPos,graph,seq,seqStr,genezilla);
  double secondScore=
    scoreSignal(rightSignal,rightPos,graph,seq,seqStr,genezilla);

  // GET HEXAMER SCORE
  double hexamerScore=scoreHexamers(begin,end,seqStr,strand);
  double hmmScore=scoreHMM(begin,end,seqStr,seq,contentType,genezilla);

  // PRODUCE OUTPUT
  if(finite(hmmScore) && finite(hexamerScore) && finite(lengthProb) &&
     finite(firstScore) && finite(secondScore) && !isnan(hmmScore))
    output(contentType,lengthProb,firstScore,secondScore,hmmScore,1);
}



void Application::output(ContentType contentType,double lengthProb,
			 double firstScore,double secondScore,
			 double hmmScore,int category)
{
  switch(contentType)
    {
    case INITIAL_EXON:
      osInitial<<lengthProb<<"\t"<<firstScore<<"\t"<<secondScore<<"\t"<<
	hmmScore<<"\t"<<category<<endl;
      break;
    case INTERNAL_EXON:
      osInternal<<lengthProb<<"\t"<<firstScore<<"\t"<<secondScore<<"\t"<<
	hmmScore<<"\t"<<category<<endl;
      break;
    case FINAL_EXON:
      osFinal<<lengthProb<<"\t"<<firstScore<<"\t"<<secondScore<<"\t"<<
	hmmScore<<"\t"<<category<<endl;
      break;
    case SINGLE_EXON:
      osSingle<<lengthProb<<"\t"<<firstScore<<"\t"<<secondScore<<"\t"<<
	hmmScore<<"\t"<<category<<endl;
      break;
    case NEG_INITIAL_EXON:
      osInitial<<lengthProb<<"\t"<<secondScore<<"\t"<<firstScore<<"\t"<<
	hmmScore<<"\t"<<category<<endl;
      break;
    case NEG_INTERNAL_EXON:
      osInternal<<lengthProb<<"\t"<<secondScore<<"\t"<<firstScore<<"\t"<<
	hmmScore<<"\t"<<category<<endl;
      break;
    case NEG_FINAL_EXON:
      osFinal<<lengthProb<<"\t"<<secondScore<<"\t"<<firstScore<<"\t"<<
	hmmScore<<"\t"<<category<<endl;
      break;
    case NEG_SINGLE_EXON:
      osSingle<<lengthProb<<"\t"<<secondScore<<"\t"<<firstScore<<"\t"<<
	hmmScore<<"\t"<<category<<endl;
      break;
    }
}



double Application::scoreSignal(SignalType signalType,int pos,
				ParseGraph &graph,const Sequence &seq,
				const BOOM::String &seqStr,
				GeneZilla &genezilla)
{
  // First, see if the signal is in the parse graph (and therefore 
  // already scored)
  int index=graph.findSignal(signalType,pos);
  if(index>-1) return graph.getIthVertex(index)->contextWindowScore();

  // Otherwise, we have to invoke the appropriate signal sensor
  // to score this putative signal
  return genezilla.getSignalSensor(signalType).getLogP(seq,seqStr,pos);
}



void Application::exonCoordsToSigPos(ContentType contentType,
				     int begin,int end,
				     int &pos1,int &pos2)
{
  switch(contentType)
    {
    case INITIAL_EXON:             // ATG - GT
      pos1=begin; pos2=end; break;
    case INTERNAL_EXON:            // AG - GT
      pos1=begin-2; pos2=end; break;
    case FINAL_EXON:               // AG - TAG
      pos1=begin-2; pos2=end-3; break;
    case SINGLE_EXON:              // ATG - TAG
      pos1=begin; pos2=end-3; break;
    case NEG_INITIAL_EXON:         // NEG_GT - NEG_ATG
      pos1=begin-2; pos2=end-3; break;
    case NEG_INTERNAL_EXON:        // NEG_GT - NEG_AG
      pos1=begin-2; pos2=end; break;
    case NEG_FINAL_EXON:           // NEG_TAG - NEG_AG
      pos1=begin; pos2=end; break;
    case NEG_SINGLE_EXON:          // NEG_TAG - NEG_ATG
      pos1=begin; pos2=end-3; break;
    }
}



void Application::sigPosToExonCoords(ContentType contentType,
				     int &begin,int &end,
				     int pos1,int pos2)
{
  switch(contentType)
    {
    case INITIAL_EXON:             // ATG - GT
      begin=pos1; end=pos2; break;
    case INTERNAL_EXON:            // AG - GT
      begin=pos1+2; end=pos2; break;
    case FINAL_EXON:               // AG - TAG
      begin=pos1+2; end=pos2+3; break;
    case SINGLE_EXON:              // ATG - TAG
      begin=pos1; end=pos2+3; break;
    case NEG_INITIAL_EXON:         // NEG_GT - NEG_ATG
      begin=pos1+2; end=pos2+3; break;
    case NEG_INTERNAL_EXON:        // NEG_GT - NEG_AG
      begin=pos1+2; end=pos2; break;
    case NEG_FINAL_EXON:           // NEG_TAG - NEG_AG
      begin=pos1; end=pos2; break;
    case NEG_SINGLE_EXON:          // NEG_TAG - NEG_ATG
      begin=pos1; end=pos2+3; break;
    }
}



void Application::getLogRatios(const BOOM::String &backgroundFilename,
			       BOOM::Map<BOOM::String,TranscriptList*> &gff,
			       const BOOM::String &substrateFilename)
{
  // Count hexamers in background sequence
  cerr<<"\tgetting background counts..."<<endl;
  BOOM::StringMap<int> backgrounds(LARGE_PRIME), codingCounts(LARGE_PRIME);
  int backgroundSS=getBackgrounds(backgroundFilename,backgrounds);
  
  // Count hexamers in coding sequence
  cerr<<"\tcounting hexamers in coding sequence..."<<endl;
  BOOM::FastaReader fastaReader(substrateFilename);
  BOOM::String defline, seqString;
  int codingSS=0;
  while(fastaReader.nextSequence(defline,seqString))
    {
      // Parse the substrate ID out of the defline
      BOOM::String contigId, remainder;
      BOOM::FastaReader::parseDefline(defline,contigId,remainder);

      // Get all the known genes on this contig from the GFF
      TranscriptList *knownGenes=gff[contigId];
      if(!knownGenes) continue;

      // Iterate through all exons
      TranscriptList::iterator cur=knownGenes->begin(), 
	end=knownGenes->end();
      for(; cur!=end ; ++cur)
	{
	  BOOM::GffTranscript &transcript=**cur;
	  BOOM::Vector<BOOM::GffExon*>::iterator cur=transcript.getExons(),
	    end=transcript.getExonsEnd();
	  for(; cur!=end ; ++cur)
	    {
	      BOOM::GffExon &exon=**cur;
	      int begin=exon.getBegin(), end=exon.getEnd();
	      int len=end-begin;
	      codingSS+=(len-6+1);
	      countHexamers(begin,end,seqString,exon.getStrand(),
			    codingCounts);
	    }
	}
    }
  double reciprocalSS=backgroundSS/double(codingSS);

  // Compute ratios
  cerr<<"\tcomputing ratios..."<<endl;
  NthOrderStringIterator iter(6,DnaAlphabet::global);
  while(!iter.done())
    {
      BOOM::String hexamer=iter.getNextString();
      const char *pHex=hexamer.c_str();
      int background=1, coding=1;
      if(backgrounds.isDefined(pHex,6))
	background+=backgrounds.lookup(pHex,6);
      if(codingCounts.isDefined(pHex,6))
	coding+=codingCounts.lookup(pHex,6);
      if(background>1 || coding>1)
	{
	  double ratio=log(coding/double(background)*reciprocalSS);
	  hexamerLogRatios.lookup(pHex,6)=ratio;
	}
    }
}



int Application::getBackgrounds(const BOOM::String &filename,
				BOOM::StringMap<int> &backgrounds)
{
  // Iterate through all sequences in the multi-fasta file:
  BOOM::FastaReader reader(filename);
  BOOM::String defline, contig;
  float gcContent;
  int sampleSize=0;
  while(reader.nextSequence(defline,contig))
    {
      // Iterate through all hexamers in this sequence & count 
      // their occurrences
      const char *p=contig.c_str();
      int len=contig.length();
      sampleSize+=(len-6+1);
      countHexamers(p,len,backgrounds);
    }

  // Return sample size
  return sampleSize;
}



double Application::scoreHexamers(int begin,int end,
				  const BOOM::String &seqStr,
				  Strand strand)
{
  // See if we need to reverse complement
  int len=end-begin;
  BOOM::String exonSeq=seqStr.substring(begin,len);
  if(strand==REVERSE_STRAND) 
    exonSeq=BOOM::ProteinTrans::reverseComplement(exonSeq);

  // Iterate through all hexamers that occur in the string,
  // accumulating their log ratios:
  const char *p=exonSeq.c_str();
  const char *endP=p+len-6+1;
  double score=0.0;
  for(; p<endP ; ++p)
    if(hexamerLogRatios.isDefined(p,6))
      score+=hexamerLogRatios.lookup(p,6);

  // Return the sum
  return score;
}



void Application::countHexamers(int begin,int end,
				const BOOM::String &seqStr,
				Strand strand,
				BOOM::StringMap<int> &counts)
{
  // May need to reverse complement:
  int len=end-begin;
  BOOM::String exonSeq=seqStr.substring(begin,len);
  if(strand==REVERSE_STRAND) 
    exonSeq=BOOM::ProteinTrans::reverseComplement(exonSeq);

  // Iterate through all hexamers:
  const char *p=exonSeq.c_str();
  const char *endP=p+len-6+1;
  for(; p<endP ; ++p)
    {
      if(!counts.isDefined(p,6)) counts.lookup(p,6)=0;

      // Increment the count for this hexamer
      ++counts.lookup(p,6);
    }
}



void Application::countHexamers(const char *p,int len,
				BOOM::StringMap<int> &counts)
{
  const char *endP=p+len-6+1;
  for(; p<endP ; ++p)
    {
      if(!counts.isDefined(p,6)) counts.lookup(p,6)=0;
      ++counts.lookup(p,6);
    }
}



int Application::countExons(TranscriptList &transcripts)
{
  // Iterate through all transcripts
  TranscriptList::iterator tCur=transcripts.begin(), 
    tEnd=transcripts.end();
  int count=0;
  for(; tCur!=tEnd ; ++tCur)
    {
      BOOM::GffTranscript &transcript=**tCur;

      // Add the number of exons contained in this transcript
      count+=transcript.getNumExons();
    }

  // Return the total number of exons
  return count;
}



void Application::getIntervals(TranscriptList &transcripts,
			       IntervalSet &intervals)
{
  // Iterate through all transcripts
  TranscriptList::iterator cur=transcripts.begin(), 
    end=transcripts.end();
  for(; cur!=end ; ++cur)
    {
      BOOM::GffTranscript &transcript=**cur;

      // Iterate through exon of this transcript
      BOOM::Vector<BOOM::GffExon*>::iterator cur=transcript.getExons(),
	end=transcript.getExonsEnd();
      for(; cur!=end ; ++cur)
	{
	  BOOM::GffExon &exon=**cur;

	  // Add this exon's interval to the list
	  intervals.insert(Interval(exon.getBegin(),exon.getEnd()));
	}
    }
}



Signal &randomExonStart(ParseGraph &orfGraph,int numVertices)
{
  // Pick a random starting position
  int start=RandomNumber(numVertices);

  // Iterate through the list until a signal of the right type 
  // is encountered
  Signal *signal;
  bool found=false;
  for(int i=0 ; i<numVertices ; ++i)
    {
      signal=orfGraph.getIthVertex((start+i)%numVertices);
      if(beginsCoding(signal->getSignalType())) {found=true; break;}
    }

  // None found?
  if(!found) throw "No coding signal found in randomExonStart()";

  // Return the signal
  return *signal;
}



void Application::outputNegatives(TranscriptList &knownGenes,
				  ostream &os,ParseGraph &orfGraph,
				  const BOOM::String &seqStr,
				  const Sequence &seq,GeneZilla &genezilla)
{
  // Cordon off the exact intervals of known exons, so we don't
  // accidentally chose any of those:
  IntervalSet usedIntervals;
  getIntervals(knownGenes,usedIntervals);

  // Generate a number of negative examples equal to the number of
  // positive examples:
  int numVertices=orfGraph.numVertices();
  int numExons=countExons(knownGenes), iter=0, maxIterations=1000000;
  for(int i=0 ; i<numExons && iter<maxIterations ; ++i, ++iter)
    {
      // Choose a random starting signal (on either strand)
      Signal &beginSignal=randomExonStart(orfGraph,numVertices);
      Strand strand=beginSignal.getStrand();
      BOOM::Set<Edge*> &edges=beginSignal.getEdgesOut();
      
      // Select a random edge to a successor signal:
      int numEdges=edges.size();
      if(numEdges==0) {--i; continue;}
      int index=RandomNumber(numEdges);
      BOOM::Set<Edge*>::iterator cur=edges.begin();
      for(int j=0 ; j<index ; ++j) ++cur;
      Edge &edge=**cur;

      // Compute the exact exon coordinates
      ContentType contentType=edge.getContentType();
      Signal &endSignal=*edge.getRight();
      int begin, end; 
      int pos1=beginSignal.getConsensusPosition();
      int pos2=endSignal.getConsensusPosition();
      sigPosToExonCoords(contentType,begin,end,pos1,pos2);
      assert(pos1<pos2);
      assert(begin<end);
      int length=end-begin;
      if(begin<50 || end+50>seqStr.length()) {--i; continue;}

      // Make sure we haven't already used this interval
      if(usedIntervals.isMember(pair<int,int>(begin,end))) {--i; continue;}

      // Compute features for this interval
      double hexamerScore=scoreHexamers(begin,end,seqStr,strand);
      double hmmScore=scoreHMM(begin,end,seqStr,seq,contentType,genezilla);
      double firstScore=beginSignal.contextWindowScore();
      double secondScore=endSignal.contextWindowScore();
      DiscreteDistribution &distribution=
	genezilla.getDistribution(contentType);
      double lengthProb=distribution.getLogP(length);
      usedIntervals.insert(pair<int,int>(begin,end));

      // Produce output
      if(finite(hmmScore) && finite(hexamerScore) && finite(lengthProb) &&
	 finite(firstScore) && finite(secondScore) && !isnan(hmmScore))
	output(contentType,lengthProb,firstScore,secondScore,hmmScore,0);
    }
  if(iter==maxIterations) 
    cout<<"WARNING: Infinite loop in outputNegatives() ... skipping"<<endl;
}



double Application::scoreHMM(int begin,int end,const BOOM::String &seqStr,
			     const Sequence &seq,ContentType contentType,
			     GeneZilla &genezilla)
{
  const int length=end-begin;
  ContentSensor &intergenicSensor=genezilla.getContentSensor(INTERGENIC);

  // Iterate through all three phases (keeping best best score of the 3)
  double bestScore=NEGATIVE_INFINITY;
  for(int phase=0 ; phase<3 ; ++phase)
    {
      ContentSensor &sensor=genezilla.getContentSensor(contentType);
      cout<<"scoring "<<contentType<<endl;
      double score=sensor.scoreSubsequence(seq,seqStr,begin,length,phase);
      cout<<"score="<<score<<endl;
      if(score>bestScore) bestScore=score;
    }
  
  double nullScore=
    intergenicSensor.scoreSubsequence(seq,seqStr,begin,length,0);
  return bestScore/nullScore;
  //return bestScore;
}
