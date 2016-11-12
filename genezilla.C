/****************************************************************
 genezilla.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "GZilla.H"
#include "genezilla.H"
#include "EdgeFactory.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Constants.H"
using namespace std;

#ifdef EXPLICIT_GRAPHS
//#error Please edit the file genezilla.H, comment out the definition of EXPLICIT_GRAPHS, issue a "make clean", and recompile this project
#endif

#ifdef FORCE_SPECIFIC_SIGNALS
#error Please edit the file genezilla.H, comment out the definition of FORCE_SPECIFIC_SIGNALS, issue a "make clean", and recompile this project
#endif

static const char *PROGRAM_NAME="genezilla";
static const char *VERSION="1.2";
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
	return -1;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
	return -1;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
	return -1;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
	return -1;
      }
    return 0;
  }



void AppMain(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"s:i:c:o:I:tP:SCDO");
  if(cmd.numArgs()!=2)
    throw BOOM::String(
    "\ngenezilla <*.iso> <*.fasta> [options]\n\
       options:\n\
          -s <N> : ignore sequences in FASTA shorter than N bases\n\
          -i <file> : load isochore predictions from file\n\
          -c <file> : load CpG island predictions from file\n\
          -o <file> : write ORF graph into file\n\
          -t : use only one left and right terminus in the ORF graph\n\
          -I <file> : dump the intergenic prefix-sum array into file\n\
          -P <dir> : use RNA-seq data from pileup\n\
          -S : omit signal scores\n\
          -C : omit content scores\n\
          -D : omit duration scores\n\
          -O : allow in-frame stop codons\n\
");
  BOOM::String isochoreFilename=cmd.arg(0);
  BOOM::String fastaFilename=cmd.arg(1);
  alphabet=DnaAlphabet::global();
  int minSeqLen=
    (cmd.option('s') ? cmd.optParm('s').asInt() : 1);
  ofstream osGraph;
  bool dumpGraph;
  if(cmd.option('o')) {
    osGraph.open(cmd.optParm('o').c_str());
    dumpGraph=true;
  }
  else dumpGraph=false;
  String psaName;
  if(cmd.option('I')) psaName=cmd.optParm('I');
  bool haveEvidence=cmd.option('P');

  // Process all the contigs
  BOOM::FastaReader fastaReader(fastaFilename);
  BOOM::String defline, seqString;
  int transcriptId=-1;
  EdgeFactory *edgeFactory=NULL;
  String evidenceDir;
  const int minSupport=1; // ###
  if(haveEvidence) {
    evidenceDir=cmd.optParm('P');
    edgeFactory=new FilteredEdgeFactory(NULL);
  }
  else edgeFactory=new EdgeFactory;
  GeneZilla genezilla(PROGRAM_NAME,VERSION,*edgeFactory,transcriptId);
#ifdef EXPLICIT_GRAPHS
  if(cmd.option('t')) genezilla.useOneTerminusOnly();
#endif
  if(cmd.option('i')) genezilla.loadIsochoreBoundaries(cmd.optParm('i'));
  if(cmd.option('c')) genezilla.loadCpGislands(cmd.optParm('c'));
  if(cmd.option('S')) genezilla.omitSignalScores();
  if(cmd.option('C')) genezilla.omitContentScores();
  if(cmd.option('D')) genezilla.omitDurationScores();
  if(cmd.option('O')) genezilla.allowPTCs();
  EvidenceFilter *evidence=NULL;
  while(fastaReader.nextSequence(defline,seqString))
    {
      Sequence seq(seqString,DnaAlphabet::global());
      if(seq.getLength()<minSeqLen) continue;

      // Parse the substrate ID out of the defline
      BOOM::String contigId, remainder;
      BOOM::FastaReader::parseDefline(defline,contigId,remainder);
      cerr<<"processing substrate "<<contigId<<"..."<<endl;

      // Load RNA evidence for this substrate
      if(haveEvidence) {
	WigBinary *wig=new WigBinary(evidenceDir+'/'+contigId+".pileup");
	RnaJunctions *junctions=
	  new RnaJunctions(evidenceDir+'/'+contigId+".junctions");
	delete evidence;
	evidence=new EvidenceFilter(minSupport,wig,junctions);
	static_cast<FilteredEdgeFactory*>(edgeFactory)->setEvidence(evidence);
	genezilla.setEvidenceFilter(evidence);
	
	genezilla.setEvidenceFilter(NULL); // ###

      }

      // Predict genes and get path
      BOOM::Stack<SignalPtr> *path=
	genezilla.processChunk(seq,seqString,isochoreFilename,contigId,
			       osGraph,dumpGraph,psaName);

      // Don't need the path; just delete it (this will also delete the
      // signal objects in it, since we use "smart pointers")
      delete path;
    }
}



