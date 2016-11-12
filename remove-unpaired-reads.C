/****************************************************************
 remove-unpaired-reads.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/Pipe.H"
using namespace std;
using namespace BOOM;


struct Read {
  static Regex idRegex;
  String lines[4];
  String id; // doesn't include the /1 or /2 at the end
  bool written;
  Read() : written(true) {}
  void load(File &);
  void write(File &);
};
Regex Read::idRegex("(.*)(\\/\d)");


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
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
  if(cmd.numArgs()!=6)
    throw String("remove-unpaired-reads <in1.fastq.gz> <in2.fastq.gz> <out1.fastq.gz> <out2.fastq.gz> <unpaired1.fastq.gz> <unpaired2.fastq.gz>");
  const String infile1=cmd.arg(0);
  const String infile2=cmd.arg(1);
  const String outfile1=cmd.arg(2);
  const String outfile2=cmd.arg(3);
  const String unpaired1=cmd.arg(4);
  const String unpaired2=cmd.arg(5);

  GunzipPipe in1(infile1), in2(infile2);
  GzipPipe out1(outfile1), out2(outfile2);
  GzipPipe singletons1(unpaired1), singletons2(unpaired2);

  Read read1, read2;
  String prevId1, prevId2;
  while(!in1.eof() && !in2.eof()) {
    read1.load(in1); read2.load(in2);
    while(read1.id!=read2.id) {
      while(read1.id<read2.id && !in1.eof()) {
	if(!read1.written) read1.write(singletons1);
	prevId1=read1.id;
	read1.load(in1);
	if(read1.id<prevId1) throw "fastq file is not sorted";
      }
      while(read2.id<read1.id && !in2.eof()) {
	if(!read2.written) read2.write(singletons2);
	prevId2=read2.id;
	read2.load(in2);
	if(read2.id<prevId2) throw "fastq file is not sorted";
      }
    }
    if(read1.id==read2.id) {
      read1.write(out1); read2.write(out2);
      prevId1=read1.id; prevId2=read2.id;
    }
  }
  if(!read1.written) read1.write(singletons1);
  if(!read2.written) read2.write(singletons2);
  while(!in1.eof()) {
    read1.load(in1);
    if(read1.id<prevId1) throw "fastq file is not sorted";
    if(!read1.written) read1.write(singletons1);
    prevId1=read1.id;
  }
  while(!in2.eof()) {
    read2.load(in2);
    if(read2.id<prevId2) throw "fastq file is not sorted";
    if(!read2.written) read1.write(singletons2);
    prevId2=read2.id;
  }

  return 0;
}



void Read::load(File &file)
{
  for(int i=0 ; i<4 ; ++i)
    lines[i]=file.getline();
  if(!idRegex.match(lines[0])) throw lines[0]+" : can't parse read id";
  id=idRegex[1];
  if(lines[0].length()>0) written=false;
  else written=true; // nothing was loaded, so nothing to write
}



void Read::write(File &file)
{
  for(int i=0 ; i<4 ; ++i) file.write(lines[i]+"\n");
  written=true;
}


