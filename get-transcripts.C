/*
 get-transcripts.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/IndexedFasta.H"
#include "BOOM/GffReader.H"
using namespace std;


class Application
{
public:
  Application();
  int main(int argc,char *argv[]);
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
    BOOM::CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=2)
      throw string("get-transcripts <*.gff> <compiled-fasta-file>");
    const BOOM::String gffFile=cmd.arg(0);
    const BOOM::String fastaFile=cmd.arg(1);

    BOOM::IndexedFasta fasta(fastaFile);
    BOOM::GffReader gffReader(gffFile);
    TranscriptList &transcripts=*gffReader.loadTranscripts();
    int numTranscripts=transcripts.size();
    for(int i=0 ; i<numTranscripts ; ++i)
      {
	BOOM::GffTranscript *transcript=transcripts[i];
	transcript->loadSequence(fasta);
	BOOM::String seq=transcript->getSequence();
	BOOM::String start=seq.substring(0,3);
	BOOM::String stop=seq.substring(seq.length()-3,3);
	if(start!="ATG" ||
	   (stop!="TAG" && stop!="TGA" && stop!="TAA"))
	  {
	    cerr<<"WARNING: skipping transcript "
		<<transcript->getTranscriptId()<<" with start codon "
		<<start<<" and stop codon "<<stop<<endl;
	    continue;
	  }
	BOOM::String transId=transcript->getTranscriptId();
	transId=transId.substitute("\"","");
	cout<<">"<<transId
	    <<" /numExons="<<transcript->getNumExons()
	    <<" /strand="<<transcript->getStrand()
	    <<" /begin="<<transcript->getBegin()
	    <<" /end="<<transcript->getEnd()
	    <<" /extent="<<(transcript->getEnd()-transcript->getBegin())
	    <<" /length="<<seq.length()
	    <<" /substrate="<<transcript->getSubstrate()
	    <<" /geneId="<<transcript->getGeneId()
	    <<" /startCodon=0"
	    <<endl;
	int len=seq.length();
	int numLines=len/60;
	int residue=len%60;
	if(residue>0) ++numLines;
	for(int i=0 ; i<numLines ; ++i)
	  {
	    int begin=i*60;
	    int end=begin+60;
	    if(end>len) end=len;
	    BOOM::String line=seq.substring(begin,end-begin);
	    cout<<line<<endl;
	  }
      }

    return 0;
  }

