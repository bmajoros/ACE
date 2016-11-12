/****************************************************************
 relative-entropy.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/Constants.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Sequence.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/File.H"
#include "ContentSensor.H"
using namespace std;
using namespace BOOM;

PureDnaAlphabet alphabet;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  void write(Array1D<double> &,File &);
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
    throw String("relative-entropy <*.model> <*.model>");
  String model1Filename=cmd.arg(0);
  String model2Filename=cmd.arg(1);

  ContentSensor *sensor1=ContentSensor::load(model1Filename);
  ContentSensor *sensor2=ContentSensor::load(model2Filename);
  bool phased=sensor1->isPhased() || sensor2->isPhased();

  double H=0.0, H1=0.0;
  int N=1+max(sensor1->getOrder(),sensor2->getOrder());
  unsigned numStrings=pow(4,N);
  Sequence seq;
  String str;
  int numPhases=(phased ? 3 : 1);

  double sumP=0.0;
  for(unsigned i=0 ; i<numStrings ; ++i) {
    for(int phase=0 ; phase<numPhases ; ++phase) {
      seq.fromInt(i,N,alphabet);
      seq.toString(alphabet,0,N,str);
      double p=sensor1->scoreSubsequence(seq,str,0,N,phase)/log(2);
      double q=sensor2->scoreSubsequence(seq,str,0,N,phase)/log(2);
      cout<<str<<"\t"<<p<<"\t"<<q<<"\t"<<p-q<<" "<<pow(2,p)*(p-q)<<endl;
      if(isFinite(p) && isFinite(q)) H+=pow(2,p)*(p-q);
      if(isFinite(p)) { sumP+=pow(2,p); H1-=pow(2,p)*p; }
    }
  }
  cout<<"sum="<<sumP<<endl;

  /*
  double p[3], q[3];
  for(unsigned i=0 ; i<numStrings ; ++i) {
    seq.fromInt(i,N,alphabet);
    seq.toString(alphabet,0,N,str);
    if(phased) {
      sensor1->scoreSingleBase(seq,str,N-1,seq[N-1],str[N-1],p[0],p[1],p[2]);
      sensor2->scoreSingleBase(seq,str,N-1,seq[N-1],str[N-1],q[0],q[1],q[2]);
      for(int j=0 ; j<3 ; ++j) {
	if(isFinite(p[j]) && isFinite(q[j])) 
	  H+=exp(p[j])*(p[j]-q[j]);
	if(isFinite(p[j]))
	  H1-=exp(p[j])*p[j];
      }
    }
    else {
      double p=sensor1->scoreSingleBase(seq,str,N-1,seq[N-1],str[N-1]);
      double q=sensor2->scoreSingleBase(seq,str,N-1,seq[N-1],str[N-1]);
      if(isFinite(p) && isFinite(q)) 
	H+=exp(p)*(p-q);
      if(isFinite(p)) 
	H1-=exp(p)*p;
    }
  }
  */
  cout<<H<<"\t"<<H1<<endl;

  return 0;
}



