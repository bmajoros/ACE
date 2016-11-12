/****************************************************************
 score-sequence.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "BOOM/FastaReader.H"
#include "BOOM/CommandLine.H"
#include "ContentSensor.H"
using namespace std;

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
  BOOM::CommandLine cmd(argc,argv,"s:i:c:");
  if(cmd.numArgs()!=2)
    throw BOOM::String("\nscore-sequence <*.model> <*.fasta>\n");
  BOOM::String modelFilename=cmd.arg(0);
  BOOM::String fastaFilename=cmd.arg(1);
  alphabet=DnaAlphabet::global();

  // Load model
  ContentSensor *model=ContentSensor::load(modelFilename);

  // Process all the sequences
  BOOM::FastaReader fastaReader(fastaFilename);
  BOOM::String defline, seqString;
  while(fastaReader.nextSequence(defline,seqString))
    {
      Sequence seq(seqString,DnaAlphabet::global());
      int length=seqString.length();
      double score=model->scoreSubsequence(seq,seqString,0,length,0);
      cout<<score<<endl;
    }
}



