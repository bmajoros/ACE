/****************************************************************
 vcf-population.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Pipe.H"
#include "BOOM/Vector.H"
#include "BOOM/Map.H"
#include "BOOM/Regex.H"
#include "BOOM/Time.H"
using namespace std;
using namespace BOOM;

struct Population {
  String name;
  Vector<int> alleleCounts;
  void increment(int allele) {
    while(alleleCounts.size()<=allele) alleleCounts.push_back(0);
    ++alleleCounts[allele];
  }
  void clearCounts() { alleleCounts.clear(); }
  int majorAllele() {
    int m=0;
    for(int i=1 ; i<alleleCounts.size() ; ++i)
      if(alleleCounts[i]>alleleCounts[m]) m=i;
    return m;
  }
  Population() {}
  Population(const String &n) : name(n) {}
};

struct Region {
  int begin, end;
  Region(int b,int e) : begin(b), end(e) {}
  bool contains(int p) { return begin<=p && p<end; }
};

struct Variant {
  String chr;
  int pos;
  String ref, id;
  Vector<String> alt;
  Variant(const String &chr,int pos,const String &ref,
	  const Vector<String> &alt,const String &id)
    : chr(chr), pos(pos), ref(ref), alt(alt), id(id) {}
  String asString() const {
    String out=id+":"+chr+":"+pos+":"+ref;
    for(Vector<String>::const_iterator cur=alt.begin(), end=alt.end() ;
	cur!=end ; ++cur) out+=String(":")+*cur;
    return out;
  }
};

struct Individual {
  String id;
  int populationIndex;
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Time timer;
  Vector<Population> populations;
  Vector<Individual> individuals;
  Map<String,int> indivToPopIndex;
  Regex gzipRegex; // *.gz
  Regex dnaRegex;
  Regex CNregex; // <CN14>
  int numIndividuals;
  void clearCounts();
  void convert(File &infile,File &outfile);
  void parseChromLine(const Vector<String> &,File &outfile);
  bool parseVariant(const Vector<String> &fields,String &chr,int &pos,
		    String &ref,Vector<String> &alt,String &id);
  bool parseVariant(const Vector<String> &,File &outfile);
  void parseVariantAndGenotypes(const Vector<String> &,const String &line,
				File &outfile);
  void processPopFile(const String &filename);
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
  : gzipRegex(".*\\.gz"), dnaRegex("^[ACGTacgt]+$"), 
    CNregex("^<CN(\\d+)>$"), numIndividuals(0)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("\nvcf-population [options] <in.vcf> <populations.txt> <out.vcf>\n\
   both input and output files can be zipped (use .gz as suffix)\n\
");
  const String infile=cmd.arg(0);
  const String popFile=cmd.arg(1);
  const String outfile=cmd.arg(2);

  // Open files
  processPopFile(popFile);
  File *vcfIn=gzipRegex.match(infile) ? 
    new Pipe(String("cat ")+infile+" | gunzip","r") : 
    new File(infile);
  File *vcfOut=gzipRegex.match(outfile) ?
    new Pipe(String("bgzip > ")+outfile,"w") :
    new File(outfile,"w");

  // Perform conversion
  timer.startCounting();
  convert(*vcfIn,*vcfOut);
  vcfIn->close(); vcfOut->close();

  return 0;
}



void Application::convert(File &infile,File &outfile)
{
  while(!infile.eof()) {
    String line=infile.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0) {
      if(fields[0]=="#CHROM") parseChromLine(fields,outfile);
      else if(fields[0][0]!='#') parseVariantAndGenotypes(fields,line,outfile);
      else outfile.print(line+"\n");
    }
    delete &fields;
  }
}



void Application::parseChromLine(const Vector<String> &fields,File &outfile)
{
  // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  200848168@1097100704
  const int numFields=fields.size();
  if(numFields<10 || fields[8]!="FORMAT") 
    throw "Error parsing #CHROM line in VCF file";
  for(int i=0 ; i<9 ; ++i) outfile.print(fields[i]+"\t");
  for(int i=0 ; i<populations.size() ; ++i) {
    outfile.print(populations[i].name);
    if(i+1<populations.size()) outfile.print("\t");
  }
  outfile.print("\n");
  int numIndiv=numFields-9;
  individuals.resize(numIndiv);
  for(int i=0 ; i<numIndiv ; ++i) {
    const String &id=fields[i+9];
    Individual &indiv=individuals[i];
    indiv.id=id;
    indiv.populationIndex=indivToPopIndex[id];
  }
  numIndividuals=numIndiv;
}



void Application::parseVariantAndGenotypes(const Vector<String> &fields,
					   const String &line,File &outfile)
{
  if(!parseVariant(fields,outfile)) return;
  const int numIndiv=fields.size()-9;
  if(numIndiv!=numIndividuals) {
    cout<<numIndiv<<"\t"<<numIndividuals<<endl;
    throw String("Number of fields in VCF line does not match the number of individuals: ")+line;
  }
  clearCounts();
  for(int i=0 ; i<numIndiv ; ++i) {
    const String &genotype=fields[i+9];
    const Individual &indiv=individuals[i];
    Population &pop=populations[indiv.populationIndex];
    Vector<String> fields; genotype.getFields(fields,"|/");
    for(Vector<String>::iterator cur=fields.begin(), end=fields.end() ;
	cur!=end ; ++cur)
      pop.increment((*cur).asInt());
  }
  for(int i=0 ; i<populations.size() ; ++i) {
    const Population &pop=populations[i];
    const int major=pop.majorAllele();
    outfile.print(major);
    if(i+1<populations.size()) outfile.print("\t");
  }
  outfile.print("\n");
}



bool Application::parseVariant(const Vector<String> &fields,
			       String &chr,int &pos,String &ref,
			       Vector<String> &alts,String &id)
{
  if(fields[6]!="PASS" && fields[6]!=".") return false;
  if(fields.size()<10) return false;
  if(fields[8]!="GT") {
    cerr<<"WARNING: GT expected in field 9 of VCF file"<<endl;
    return false;
  }
  if(fields[3].contains("<") || fields[4].contains("<")) return false;
  chr=fields[0];
  pos=fields[1].asInt()-1; // VCF files are 1-based
  id=fields[2];
  if(id==".") id=chr+"@"+pos;
  ref=fields[3];
  const String &alt=fields[4];
  if(!dnaRegex.match(ref)) return false;
  alt.getFields(alts,","); // in case it's a multi-allelic locus
  for(Vector<String>::iterator cur=alts.begin(), end=alts.end() ; cur!=end ;
      ++cur) {
    String &alt=*cur;
    if(CNregex.match(alt)) { // CNV (copy-number variant): <CN12>
      const int n=CNregex[1];
      alt="";
      for(int i=0 ; i<n ; ++i) alt+=ref;
    }
    else if(!dnaRegex.match(alt)) return false;
  }
  return true;
}



bool Application::parseVariant(const Vector<String> &fields,File &outfile)
{
  String chr, ref, id;
  Vector<String> alt;
  int pos;
  if(!parseVariant(fields,chr,pos,ref,alt,id)) return false;
  for(int i=0 ; i<9 ; ++i) outfile.print(fields[i]+"\t");
  return true;
}



void Application::processPopFile(const String &filename)
{
  Map<String,int> popMap;
  File file(filename);
  while(!file.eof()) {
    String line=file.getline();
    line.trimWhitespace(); if(line.empty()) continue;
    Vector<String> fields; line.getFields(fields);
    if(fields.size()<2) continue;
    String indiv=fields[0], pop=fields[1];
    if(!popMap.isDefined(pop)) {
      popMap[pop]=populations.size();
      populations.push_back(Population(pop));
    }
    indivToPopIndex[indiv]=popMap[pop];
  }
}



void Application::clearCounts()
{
  for(Vector<Population>::iterator cur=populations.begin(), 
	end=populations.end() ; cur!=end ; ++cur)
    (*cur).clearCounts();
}



