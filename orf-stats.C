/****************************************************************
 orf-stats
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
#include "BOOM/Array1D.H"
#include "NthOrderStringIterator.H"
#include <assert.h>
using namespace std;



/****************************************************************
                         defines and pragmas
 */
#ifdef EXPLICIT_GRAPHS
#error Please edit the file genezilla.H, comment out the line defining EXPLICIT_GRAPHS, issue a "make clean" and then recompile this project
#endif
const int LARGE_PRIME=3121;
int minOrfLen;



/****************************************************************
                               globals
 */
static const char *PROGRAM_NAME="get-classifier-examples";
static const char *VERSION="1.0";
Alphabet alphabet;
int frame; // ### CAUTION: this is required by older code; to be removed
GarbageCollector garbageCollector;
BOOM::String contigId;
bool coordsOnly;
bool applyCutoff;

//------------------------------------------------------------
class Application
{
public:
  Application();
  int main(int argc,char *argv[]);
};

//------------------------------------------------------------
struct ORF
{
  int begin, end, phase;
  double score;
  Strand strand;
  ContentType contentType;
  ORF(int b,int e,int p,double sc,Strand s,ContentType t)
    : begin(b), end(e), phase(p), score(sc), strand(s), contentType(t) {}
  bool overlaps(const ORF &o) {return begin<o.end && end>o.begin;}
  bool contains(const ORF &o) {return begin<=o.begin && end>=o.end;}
  int length() {return end-begin;}
};

//------------------------------------------------------------
class OrfScan : public GeneZilla
{
  BOOM::Array1D<ORF*> orfBases;
  BOOM::Set<ORF*> orfs;

  void addOrf(ORF *);
  virtual void updateAccumulators(const Sequence &seq,
				  const BOOM::String &str,
				  int pos,Symbol base,char c);
  virtual void selectCodingPred(int newConsPos,SignalQueue &queue,
			       Strand strand,double bestScore[3],
			       SignalPtr bestPred[3],ContentType,
			       SignalType toType,SignalPtr);
  virtual BOOM::Stack<SignalPtr> * mainAlgorithm(const Sequence &,
					       const BOOM::String &);
public:
  virtual BOOM::Stack<SignalPtr> * mainAlgorithm(const Sequence &,
					       const BOOM::String &,
					       const BOOM::String &
					         isochoreFilename,
					       const BOOM::String &substrateId);
  OrfScan(const BOOM::String &PROGRAM_NAME,const BOOM::String &VERSION,
	   EdgeFactory &,int &transcriptId);  
};

int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
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
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    BOOM::CommandLine cmd(argc,argv,"cT");
    if(cmd.numArgs()!=5)
      throw string("\n\
orf-stats [-ct] <*.iso> <*.fasta> <min-substrate-len> <min-orf-len> <max-total-length>\n\
    -c : output (c)oordinates only\n\
    -T : ignore signal sensor (T)hresholds\n\
\n");
    BOOM::String isoFilename=cmd.arg(0);
    BOOM::String fastaFilename=cmd.arg(1);
    int minLen=cmd.arg(2).asInt();
    minOrfLen=cmd.arg(3).asInt();
    int maxLen=cmd.arg(4).asInt();
    coordsOnly=cmd.option('c');
    applyCutoff=cmd.option('t');
    alphabet=DnaAlphabet::global;

    // Process all the contigs
    BOOM::FastaReader fastaReader(fastaFilename);
    BOOM::String defline, seqString;
    int transcriptId=1000;
    float gcContent;
    int totalLen=0;
    while(fastaReader.nextSequence(defline,seqString))
      {
	Sequence seq(seqString,DnaAlphabet::global);
	int len=seq.getLength();
	if(len<minLen) continue;
	
	// Parse the substrate ID out of the defline
	BOOM::String remainder;
	BOOM::FastaReader::parseDefline(defline,contigId,remainder);
	cerr<<"processing contig "<<contigId<<"..."<<endl;
	
	// Run the orf-finder
	cerr<<"identifying open reading frames..."<<endl;
	EdgeFactory edgeFactory;
	OrfScan orfscan(PROGRAM_NAME,VERSION,edgeFactory,transcriptId);
	//orfscan.dontUseBayesModel();
	BOOM::Stack<SignalPtr> *path=orfscan.mainAlgorithm(seq,seqString,
							 isoFilename,contigId);
	delete path;

	totalLen+=len;
	if(!coordsOnly)
	  cout<<"TOTAL LENGTH: "<<totalLen<<endl;
	cerr<<"TOTAL LENGTH: "<<totalLen<<endl;
	if(totalLen>=maxLen) break;
      }

    return 0;
  }



OrfScan::OrfScan(const BOOM::String &PROGRAM_NAME,const BOOM::String &VERSION,
		 EdgeFactory &factory,int &transcriptId)
  : GeneZilla(PROGRAM_NAME,VERSION,factory,transcriptId),
    orfBases(0)
{
}



BOOM::Stack<SignalPtr> *OrfScan::mainAlgorithm(const Sequence &seq,
					       const BOOM::String &str)
{
  const char *charPtr=str.c_str();
  instantiateLeftTermini();
  intergenicSums.resize(seqLen);
  orfBases.resize(seqLen);
  orfBases.setAllTo(NULL);
  computeIntergenicSums(seq,str,charPtr);
  for(int pos=0 ; pos<seqLen ; ++pos, ++charPtr)
    {
      Symbol base=seq[pos];
      BOOM::Vector<SignalSensor*>::iterator cur=signalSensors.begin(),
	end=signalSensors.end();
      for(; cur!=end ; ++cur )
	{
	  SignalSensor &sensor=**cur;
	  if(pos+sensor.getContextWindowLength()>seqLen) continue;

	  SignalPtr signal=(applyCutoff ?
			    sensor.detect(seq,str,pos) :
			    sensor.detectWithNoCutoff(seq,str,pos));
	  if(signal)
	    {
	      if(!coordsOnly)
		cout<<"signal "<<signal->getSignalType()
		    <<" "<< signal->contextWindowScore()<<endl;
	      linkBack(str,signal);
	      enqueue(signal);
	    }
	}
      if(pos>1) handleStopCodons(str,pos-2);
      updateAccumulators(seq,str,pos,base,*charPtr);
    }
  double parseScore;
  BOOM::Stack<SignalPtr> *path=instantiateRightTermini(str,seqLen,parseScore);

  BOOM::Set<ORF*>::iterator cur=orfs.begin(), end=orfs.end();
  for(; cur!=end ; ++cur) 
    {
      ORF *orf=*cur;
      cout<<contigId<<"\torf-stats\t"
	  <<contentTypeNiceString(orf->contentType)<<"\t"
	  <<orf->begin<<'\t'
	  <<orf->end<<'\t'<<orf->score<<'\t'<<orf->strand<<'\t'
	  <<orf->phase<<endl;    
    }

  return path;
}



void OrfScan::selectCodingPred(int newConsPos,SignalQueue &queue,
			       Strand strand,double bestScore[3],
			       SignalPtr bestPred[3],
			       ContentType contentType,
			       SignalType toType,SignalPtr signal)
{
  BOOM::Iterator<SignalPtr> &cur=queue.begin(), &end=queue.end();
  for(; cur!=end ; ++cur)
    {
      SignalPtr pred=*cur;
      Propagator &predProp=pred->getPropagator(contentType);
      SignalType predType=pred->getSignalType();
      int oldPos=pred->posOfBaseFollowingConsensus();
      int length=newConsPos-oldPos;
      double lengthScore=
	queue.getDistribution().getLogP(length+pred->getConsensusLength());
      int frameDelta=length % 3;
      for(int oldPhase=0 ; oldPhase<3 ; ++oldPhase)
	{
	  if(isinf(predProp[oldPhase])) continue;
	  int newPhase=(strand==FORWARD_STRAND ?
			(oldPhase+frameDelta) % 3 :
			posmod(oldPhase-frameDelta));
	  if(recentlyEclipsedPhases[newPhase]) continue;
	  double transScore=
	    scoreIntronPhases(predType,toType,oldPhase,newPhase);
	  double predScore=predProp[oldPhase];
	  double &bestPredScore=bestScore[newPhase];
	  double numerator=predScore-pred->posteriorInductiveScore(oldPhase);
	  double denominator=
	    intergenicSums[signal->getContextWindowPosition()-1]-
	    intergenicSums[pred->getContextWindowPosition()+
			   pred->getContextWindowLength()-1];
	  double exonScore=numerator-denominator;
	  if(coordsOnly)
	    {
	      if(length>=minOrfLen)
		addOrf(new ORF(oldPos+1,newConsPos,oldPhase,exonScore,strand,
			       contentType));
	      /*cout<<contigId<<'\t'<<"\tlong\tORF\t"<<(oldPos+1)<<'\t'
		<<newConsPos<<'\t'<<exonScore<<'\t'<<strand<<'\t'
		<<oldPhase<<endl;*/
	    }
	  else
	    cout<<"content "<<contentType<<" "<<predScore<<" "<<length<<endl;
	  if(finite(predScore) && predScore>bestPredScore)
	    {
	      bestPredScore=predScore;
	      bestPred[newPhase]=pred;
	    }
	}
    }
}



BOOM::Stack<SignalPtr> *OrfScan::mainAlgorithm(const Sequence &seq,
					     const BOOM::String &str,
					     const BOOM::String &
					        isochoreFilename,
					     const BOOM::String &substrateId)
{
  return GeneZilla::mainAlgorithm(seq,str,isochoreFilename,substrateId);
}



void OrfScan::updateAccumulators(const Sequence &seq,
				 const BOOM::String &str,
				 int pos,Symbol base,char c)
{
  BOOM::Vector<SignalQueue*>::iterator cur=signalQueues.begin(),
    end=signalQueues.end();
  for(; cur!=end ; ++cur)
    {
      SignalQueue &queue=**cur;
      ContentSensor &contentSensor=queue.getContentSensor();
      if(contentSensor.isCoding())
	{
	  double scorePhase0, scorePhase1, scorePhase2;
	  contentSensor.scoreSingleBase(seq,str,pos,base,c,scorePhase0,
					scorePhase1,scorePhase2);
	  queue.addToAccumulator(scorePhase0,scorePhase1,scorePhase2,pos);
	}
      else
	{
	  //### Hmmmm... I think the INTERGENIC model should be turned off here...
	  double score=contentSensor.scoreSingleBase(seq,str,pos,base,c);
	  queue.addToAccumulator(score);
	}
    }
  /*
  BOOM::Vector<SignalQueue*>::iterator cur=signalQueues.begin(),
    end=signalQueues.end();
  for(; cur!=end ; ++cur)
    {
      SignalQueue &queue=**cur;
      ContentSensor &contentSensor=queue.getContentSensor();
      if(contentSensor.isCoding())
	queue.addToAccumulator(0,0,0,pos);
      else
	queue.addToAccumulator(0);
    }
  */
}



void OrfScan::addOrf(ORF *newOrf)
{
  int begin=newOrf->begin, end=newOrf->end, length=newOrf->length();
  for(int i=begin ; i<end ; ++i)
    {
      ORF *otherOrf=orfBases[i];
      if(otherOrf && otherOrf->length()>length)
	{
	  delete newOrf;
	  return;
	}
    }
  for(int i=begin ; i<end ; ++i)
    {
      ORF *otherOrf=orfBases[i];
      if(otherOrf)
	{
	  int begin=otherOrf->begin, end=otherOrf->end;
	  for(int i=begin ; i<end ; ++i) orfBases[i]=NULL;
	  orfs.remove(otherOrf);
	  delete otherOrf;
	}
      orfBases[i]=newOrf;
    }
  orfs.insert(newOrf);
}


