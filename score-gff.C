/****************************************************************
 score-gff
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "GeneZilla.H"
#include "genezilla.H"
#include "EdgeFactory.H"
#include "BOOM/FastaReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/Constants.H"
using namespace std;

#ifndef EXPLICIT_GRAPHS
#error Please edit the file genezilla.H, uncomment the definition of EXPLICIT_GRAPHS, issue a "make clean", and recompile this project
#endif

#ifndef FORCE_SPECIFIC_SIGNALS
#error Please edit the file genezilla.H, uncomment the definition of FORCE_SPECIFIC_SIGNALS, issue a "make clean", and recompile this project
#endif

static const char *PROGRAM_NAME="score-gff";
static const char *VERSION="1.0";
Alphabet alphabet;
int frame; // ### CAUTION: this is required by older code; to be removed

void AppMain(int argc,char *argv[]);

int main(int argc,char *argv[])
  {
    try
      {
	AppMain(argc,argv);
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



void AppMain(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"s:r:");
  if(cmd.numArgs()!=3)
    throw BOOM::String("\nscore-gff <*.iso> <*.fasta> <*.gff>\n");
  BOOM::String isochoreFilename=cmd.arg(0);
  BOOM::String fastaFilename=cmd.arg(1);
  BOOM::String gffFilename=cmd.arg(2);
  alphabet=DnaAlphabet::global;
  
  // Load GFF
  BOOM::GffReader gffReader(gffFilename);
  BOOM::Vector<BOOM::GffTranscript*> &transcripts=*gffReader.loadTranscripts();

  // Construct parse graph and extract GFF path
  EdgeFactory edgeFactory;
  int transcriptId=999;
  GeneZilla genezilla(PROGRAM_NAME,VERSION,edgeFactory,transcriptId);
  GffPathFromParseGraph gpfpg(genezilla);
  BOOM::Vector<int> *signalCoords=gpfpg.getSignalCoordinates(transcripts);
  genezilla.forceSignalCoords(*signalCoords);
  float gcContent;
  Sequence *seq;
  BOOM::String *seqStr;
  genezilla.parse(fastaFilename,isochoreFilename,seq,seqStr,gcContent);
  BOOM::Vector<BOOL> found;
  BOOM::Vector<SignalPtr> &path=*
    genezilla.getPathFromGff(transcripts,*seq,*seqStr,found);

  // Print out scores
  BOOM::String substrateId=genezilla.getSubstrateId();
  int n=path.size();
  bool allFound=true;
  cout<<"# SIGNALS:"<<endl;
  for(int i=0 ; i<n ; ++i)
    {
      SignalPtr signal=path[i];
      BOOM::String signalName=signalTypeToName(signal->getSignalType());
      int begin=signal->getConsensusPosition();
      int end=begin+signal->getConsensusLength();
      Strand strand=signal->getStrand();
      float normScore=exp(signal->contextWindowScore()/
			  signal->getContextWindowLength());
      BOOM::String foundMsg=found[i] ? "found" : "MISSING";
      cout<<substrateId<<"\t"<<signalName
	  <<"\t"<<foundMsg<<"\t"<<begin<<"\t"<<end<<"\t"
	  <<normScore<<"\t"<<strand
	  <<"\t.\tlogp="
	  <<signal->contextWindowScore()
	  <<";"<<endl;
      if(!found[i]) allFound=false;
    }
  if(allFound)
    {
      cout<<"# EXONS:"<<endl;
      int phase=-1;
      for(int i=0 ; i<n-1 ; ++i)
	{
	  SignalPtr thisSignal=path[i], nextSignal=path[i+1];
	  SignalType thisType=thisSignal->getSignalType(), 
	    nextType=nextSignal->getSignalType();
	  Strand strand=thisSignal->getStrand();
	  if(thisType==ATG) {phase=0;}
	  else if(thisType==NEG_TAG) {phase=2;}
	  if(phase<0) throw "GFF must consist of full gene models";
	  
	  if(beginsCoding(thisType))
	    {
	      ContentType contentType;
	      int begin=thisSignal->posOfBaseFollowingConsensus();
	      int end=nextSignal->getConsensusPosition();
	      if(thisType==ATG || thisType==NEG_TAG) begin-=3;
	      if(nextType==TAG || nextType==NEG_ATG) end+=3;
	      int frameDelta=(end-begin)%3;
	      int newPhase=(strand==FORWARD_STRAND ?
			    (phase+frameDelta) % 3 :
			    posmod(phase-frameDelta));
	      double score=
		genezilla.scoreExon(thisSignal,nextSignal,phase,contentType);
	      cout<<substrateId<<"\t"
		  <<contentTypeNiceString(contentType)<<"\tfound\t"
		  <<begin<<"\t"<<end<<"\t"
		  <<score<<"\t"<<strand<<"\t"<<phase<<endl;
	      phase=newPhase;
	    }
	}
    }
  else cout<<"# SOME SIGNALS WERE NOT FOUND -- CANNOT EVALUATE GENE MODEL"<<endl;
}
