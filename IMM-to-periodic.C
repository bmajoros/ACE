/****************************************************************
 IMM-to-periodic.C
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"
#include "ThreePeriodicIMM.H"
using namespace std;
using namespace BOOM;

DnaAlphabet alphabet;

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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("IMM-to-periodic <in.IMM> <out.3PIMM>");
  const String infile=cmd.arg(0);
  const String outfile=cmd.arg(1);

  IMM *imm=new IMM(infile);
  const ContentType contentType=imm->getContentType();
  ThreePeriodicIMM *imm3p=new ThreePeriodicIMM(FORWARD_STRAND,contentType);
  for(int i=0 ; i<3 ; ++i) imm3p->getChain(i)=imm;
  //  imm3p->reverseComplement();
  imm3p->save(outfile);

  return 0;
}

