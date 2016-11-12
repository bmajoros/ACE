/****************************************************************
 entropy.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/SequenceEntropy.H"
#include "BOOM/FastaReader.H"
#include "BOOM/SummaryStats.H"
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
    if(cmd.numArgs()!=1)
      throw string("entropy <*.fasta>");
    BOOM::String filename=cmd.arg(0);

    // Process all sequences in FASTA file
    BOOM::Vector< BOOM::Vector<double> > samples(9);
    BOOM::FastaReader reader(filename);
    BOOM::String def, seq;
    double maxEntropy;
    while(reader.nextSequence(def,seq))
      {
	if(seq.length()<40) continue;
	for(int order=0 ; order<=8 ; ++order)
	  {
	    double H=BOOM::SequenceEntropy::entropy(seq,order,maxEntropy);
	    if(H==0) continue;
	    H/=maxEntropy;
	    samples[order].push_back(H);
	  }
      }

    // Compute means
    BOOM::Vector<double> curves[3];
    for(int order=0 ; order<=8 ; ++order)
      {
	BOOM::SummaryStats stats(samples[order]);
	double sigma=stats.getStdDev();
	double mu=stats.getMean();
	curves[0].push_back(mu-sigma);
	curves[1].push_back(mu);
	curves[2].push_back(mu+sigma);
      }

    // Produce output
    for(int i=0 ; i<3 ; ++i)
      {
	BOOM::Vector<double> &v=curves[i];
	for(int order=0 ; order<=8 ; ++order)
	  cout<<order<<"\t"<<v[order]<<endl;
	cout<<endl;
      }

    return 0;
  }

