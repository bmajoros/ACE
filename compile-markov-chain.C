/********************************************************************
 compile-markov-chain.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.

 Compiles a Markov chain or IMM so that it uses a transition matrix
 rather than a hash table. 

 Results on a 1.8Mb sequence:
    using only compiled content sensors = 35 sec
    using hash-table content sensors    = 1 min 15 sec
    gene predictions were the same

    ==> compiling yields about a two-fold improvement in speed

********************************************************************/

#include <string>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/Regex.H"
#include "MarkovChain.H"
#include "FastMarkovChain.H"
#include "MarkovChainCompiler.H"

Alphabet alphabet;

class Application
{
  BOOM::Regex filenameRegex;
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
  : filenameRegex("([^/]+).model$")
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    BOOM::CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=1)
      throw string("\n\
compile-markov-chain <*.model>   \n\
\n\
NOTES:\n\
       (1) output is written to *.binmod\n\
       (2) input file can be an MC, 3PMC, IMM, or 3PIMM\n\
");
    BOOM::String infile=cmd.arg(0);
    if(!filenameRegex.search(infile))
      throw BOOM::String("Can't parse filename: ")+infile;
    BOOM::String outfile=filenameRegex[1]+".binmod";

    // Misc. initialization
    alphabet=DnaAlphabet::global();
    ContentSensor *model=ContentSensor::load(infile);

    // Perform the compilation
    ContentSensor *fmc=model->compile();

    // Save the resulting compiled model
    BOOM::File file(outfile,"w");
    fmc->save(file);

    return 0;
  }

