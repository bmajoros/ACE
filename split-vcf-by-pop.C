/****************************************************************
 split-vcf-by-pop.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Vector.H"
#include "BOOM/Map.H"
#include "BOOM/Set.H"
#include "BOOM/Regex.H"
using namespace std;
using namespace BOOM;

struct Individual {
  String id;
  File *file;
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Regex dnaRegex;
  bool prependChr;
  bool SNPsOnly;
  bool variableOnly;
  bool quiet;
  Vector<Individual> individuals;
  Set<String> populations;
  Map<String,String> indivToPop;
  Vector<File*> files;
  Map<String,File*> popToFile;
  void loadPops(const String &filename);

  void filter();
  void parseChromLine(const Vector<String> &);
  void parseVariant(const Vector<String> &);
  void parseGenotype(const String &,int gt[2]);
  bool keep(const String &chr,int pos);
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
  catch(...)
    {cerr << "Unknown exception caught in main" << endl;}
  return -1;
}



Application::Application()
  : dnaRegex("[^ACGTacgt]")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"cqsv");
  if(cmd.numArgs()!=2)
    throw String("\n\
cat in.vcf.fz | gunzip | split-vcf-by-pop <pop.txt> <prefix>\n\
   pop.txt has 2 columns: individual ID, and population name\n\
   output files will be named: <prefix>-<pop>.vcf\n\
   -c : prepend \"chr\" before chromosome names\n\
   -s : SNPs only\n\
   -v : variable sites only\n\
   -q : quiet (no warnings)\n\
");
  const String popFile=cmd.arg(0);
  const String prefix=cmd.arg(1);
  variableOnly=cmd.option('v');
  prependChr=cmd.option('c');
  SNPsOnly=cmd.option('s');
  quiet=cmd.option('q');

  // Read population assignments
  loadPops(popFile);

  // Open VCF files
  for(Set<String>::iterator cur=populations.begin(), end=populations.end() ;
      cur!=end ; ++cur) {
    String pop=*cur;
    String filename=prefix+"-"+pop+".vcf";
    File *f=new File(filename,"w");
    popToFile[pop]=f;
    files.push_back(f);
  }

  // Perform conversion
  filter();

  // Clean up
  for(Vector<File*>::iterator cur=files.begin(), end=files.end() ; 
      cur!=end ; ++cur) (*cur)->close();

  return 0;
}



void Application::loadPops(const String &filename)
{
  File f(filename);
  while(!f.eof()) {
    String line=f.getline();
    if(line.isEmpty()) continue;
    Vector<String> fields; line.getFields(fields);
    if(fields.size()<2) continue;
    String indivID=fields[0], pop=fields[1];
    indivToPop[indivID]=pop;
    populations+=pop;
  }
}



void Application::filter()
{
  String headerLines;
  while(1) {
    String line;
    line.getline(cin);
    if(line.isEmpty()) break;
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0) {
      if(fields[0]=="#CHROM") {
	for(Vector<File*>::iterator cur=files.begin(), end=files.end() ; 
	    cur!=end ; ++cur) (*cur)->print(headerLines);
	parseChromLine(fields);
      }
      else if(fields[0][0]!='#') parseVariant(fields);
      else if(fields[0][0]=='#') headerLines+=line+"\n";
    }
    delete &fields;
  }
}



void Application::parseChromLine(const Vector<String> &fields)
{
  // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  200848168@1097100704
  const int numFields=fields.size();
  if(numFields<10 || fields[8]!="FORMAT") 
    throw "Error parsing #CHROM line in VCF file";
  String base=fields[0];
  for(int i=1 ; i<9 ; ++i) base+=String("\t")+fields[i];
  for(Vector<File*>::iterator cur=files.begin(), end=files.end() ; 
      cur!=end ; ++cur) (*cur)->print(base);
  int numIndiv=numFields-9;
  individuals.resize(numIndiv);
  for(int i=0 ; i<numIndiv ; ++i) {
    Individual &indiv=individuals[i];
    String id=indiv.id=fields[i+9];
    String pop=indivToPop[id];
    File *f=indiv.file=popToFile[pop];
    f->print(String("\t")+id);
  }
  for(Vector<File*>::iterator cur=files.begin(), end=files.end() ; 
      cur!=end ; ++cur) (*cur)->print(String("\n"));
}



void Application::parseVariant(const Vector<String> &fields)
{
  if(variableOnly) {
    bool seenZero=false, seenOne=false;
    for(int i=9 ; i<fields.size() ; ++i) {
      const String &field=fields[i];
      if(field.length()!=3) throw String("Bad field in VCF: ")+field;
      if(field[0]=='0' || field[3]=='0') seenZero=true;
      if(field[0]=='1' || field[3]=='1') seenOne=true;
    }
    if(!seenZero || !seenOne) return;
  }
  if(fields.size()<10 || fields[6]!="PASS" || fields[8]!="GT") return;
  const String chr=fields[0];
  if(prependChr) chr=String("chr")+chr;
  const int pos=fields[1].asInt();
  const String id=fields[2];
  if(id==".") id=chr+"@"+pos;
  const String ref=fields[3];
  const String alt=fields[4];
  if(dnaRegex.search(ref) || dnaRegex.search(alt)) return; // nonstandard chars
  if(ref.contains("<") || alt.contains("<")) {
    if(!quiet) cerr<<"skipping "<<id<<" : nonstandard variant"<<endl;
    return;
  }
  if(SNPsOnly && (ref.size()!=1 || alt.size()!=1)) return;
  String base=fields[0];
  for(int i=1 ; i<9 ; ++i) base+=String("\t")+fields[i];
  for(Vector<File*>::iterator cur=files.begin(), end=files.end() ; 
      cur!=end ; ++cur) (*cur)->print(base);
  const int numIndiv=fields.size()-9;
  for(int i=0 ; i<numIndiv ; ++i) {
    Individual &indiv=individuals[i];
    indiv.file->print(String("\t")+fields[i+9]);
  }
  for(Vector<File*>::iterator cur=files.begin(), end=files.end() ; 
      cur!=end ; ++cur) (*cur)->print(String("\n"));
}



void Application::parseGenotype(const String &s,int gt[2])
{
  if(s.length()!=3) throw String("can't parse genotype in vcf file: ")+s;
  if(s[1]!='|') throw "vcf file does not appear to be phased";
  gt[0]=s[0]-'0';
  gt[1]=s[2]-'0';
}





