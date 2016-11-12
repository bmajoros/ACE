/*
 WMM-add-pseudocounts.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <string>
#include <iostream>
#include <math.h>
#include "GarbageCollector.H"
#include "WMM.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"

Alphabet alphabet;

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
      throw string("\n\
WMM-add-pseudocounts <pseudocount>  <*.model>\n\
\n\
where <pseudocount> is the percentage of the lowest value in each position\n\
                    to use as the value inserted in place of -inf entries\n\
\n\
example: WMM-add-pseudocounts 0.5 start-codons0-100.model\n\
\n\
");
    float percentage=cmd.arg(0).asFloat();
    BOOM::String modelFile=cmd.arg(1);

    // Load the WMM
    alphabet=DnaAlphabet::global();
    GarbageIgnorer GC;
    WMM wmm(GC,modelFile);

    // Process each column
    Symbol N=alphabet.lookup('N');
    BOOM::Array2D<float> &matrix=wmm.getMatrix();
    int height=matrix.getSecondDim();
    int consensusOffset=wmm.getConsensusOffset();
    int consensusLength=wmm.getConsensusLength();
    int firstSkip=consensusOffset, lastSkip=consensusOffset+consensusLength;
    int windowLength=wmm.getContextWindowLength();
    float smallestValue=0;
    for(int i=0 ; i<windowLength ; ++i)
      {
	// Skip the consensus positions
	if(i>=firstSkip && i<lastSkip) continue;

	// Find smallest (finite) entry
	BOOM::Array2D<float>::RowIn2DArray<float> column=matrix[i];
	for(int j=0 ; j<height ; ++j)
	  {
	    if(j==N) continue;
	    float &entry=column[j];
	    if(!isinf(entry) && entry<smallestValue)
	      smallestValue=entry;
	  }

	// Replace -inf values with the new value
	float newEntry=percentage*smallestValue, total=0;
	for(int j=0 ; j<height ; ++j)
	  {
	    if(j==N) continue;
	    float &entry=column[j];
	    if(isinf(entry)) entry=newEntry;
	    total+=exp(entry);
	  }

	// Re-normalize matrix entries to sum to 1
	for(int j=0 ; j<height ; ++j) 
	  {
	    if(j==N) continue;
	    column[j]=log(exp(column[j])/total);
	  }
	
	/*// Verify that they do now sum to 1
	total=0;
	for(int j=0 ; j<height ; ++j)
	  {
	    if(j==N) continue;
	    float &entry=column[j];
	    total+=exp(entry);
	  }
	  cout<<i<<" total: "<<total<<endl;*/
      }
    wmm.save(modelFile);

    return 0;
  }

