/****************************************************************
 vcf-to-tvf.C
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
  Vector<String> genotypes;
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Time timer;
  Map<String,Vector<Region> > regions;
  Vector<Individual> individuals;
  Vector<Variant> variants;
  Regex gzipRegex; // *.gz
  Regex dnaRegex;
  Regex CNregex; // <CN14>
  bool wantFilter;
  bool prependChr;
  bool SNPsOnly;
  bool variableOnly;
  bool quiet;
  bool filterByIndiv; // user wants only specific individuals
  Set<String> keepIndiv;
  Set<String> males;
  int numIndividuals;
  bool knowMales;
  void loadIndivList(const String &filename);
  void loadRegions(const String &filename);
  void convert(File &infile,File &outfile);
  void preprocess(File &infile);
  void parseChromLine(const Vector<String> &);
  bool parseVariant(const Vector<String> &fields,String &chr,int &pos,
		    String &ref,Vector<String> &alt,String &id);
  bool parseVariant(const Vector<String> &);
  void parseVariantAndGenotypes(const Vector<String> &,const String &line);
  bool keep(const String &chr,int pos);
  void output(File &outfile);
  bool variableSite(const Vector<String> &fields);
  void loadGender(const String &filename);
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
  CommandLine cmd(argc,argv,"ci:f:qsvm:y:");
  if(cmd.numArgs()!=2)
    throw String("\nvcf-to-tvf [options] <in.vcf> <out.tvf>\n\
   both input and output files can be zipped (use .gz as suffix)\n\
   -f regions.bed : keep only variants in these regions\n\
        (coordinates are 0-based, end is not inclusive)\n\
   -c : prepend \"chr\" before chromosome names\n\
   -i <file> : keep only these individuals\n\
   -s : SNPs only\n\
   -v : variable sites only\n\
   -q : quiet (no warnings)\n\
   -m <tempfile> : small memory footprint (may be slow)\n\
   -y <gender.txt> : when not all individuals are given, assume females\n\
                     are missing (as for the Y chromosome)\n\
                     format: name followed by gender (male/female)\n\
");
  const String infile=cmd.arg(0);
  const String outfile=cmd.arg(1);
  variableOnly=cmd.option('v');
  wantFilter=cmd.option('f');
  prependChr=cmd.option('c');
  SNPsOnly=cmd.option('s');
  quiet=cmd.option('q');
  if(cmd.option('y')) loadGender(cmd.optParm('y'));
  const bool smallmem=cmd.option('m');
  if(smallmem) throw "option -m is not currently supported";
  filterByIndiv=cmd.option('i');
  if(filterByIndiv) loadIndivList(cmd.optParm('i'));
  if(!cmd.option('y') && !cmd.option('i')) throw "-i or -y must be specified";

  // Load regions to filter by
  if(wantFilter) loadRegions(cmd.optParm('f'));

  // Open files
  File *vcf=gzipRegex.match(infile) ? 
    new Pipe(String("cat ")+infile+" | gunzip","r") : 
    new File(infile);
  File *tvf=gzipRegex.match(outfile) ?
    new Pipe(String("bgzip > ")+outfile,"w") :
    new File(outfile,"w");

  // Perform conversion
  timer.startCounting();
  convert(*vcf,*tvf);
  vcf->close(); tvf->close();

  return 0;
}



void Application::loadRegions(const String &filename)
{
  File file(filename);
  while(!file.eof()) {
    String line=file.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0) {
      if(fields.size()<3) throw filename+" - can't parse bed file: "+line;
      const String &chr=fields[0];
      const int begin=fields[1].asInt();
      const int end=fields[2].asInt();
      if(!regions.isDefined(chr)) regions[chr]=Vector<Region>();
      regions[chr].push_back(Region(begin,end));
    }
    delete &fields;
  }
}



void Application::convert(File &infile,File &outfile)
{
  while(!infile.eof()) {
    String line=infile.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0) {
      if(fields[0]=="#CHROM") parseChromLine(fields);
      else if(fields[0][0]!='#') parseVariantAndGenotypes(fields,line);
    }
    delete &fields;
  }
  output(outfile);
}



void Application::parseChromLine(const Vector<String> &fields)
{
  // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  200848168@1097100704
  const int numFields=fields.size();
  if(numFields<10 || fields[8]!="FORMAT") 
    throw "Error parsing #CHROM line in VCF file";
  int numIndiv=numFields-9;
  individuals.resize(numIndiv);
  for(int i=0 ; i<numIndiv ; ++i)
    individuals[i].id=fields[i+9];
}



void Application::parseVariantAndGenotypes(const Vector<String> &fields,
					   const String &line)
{
  if(!parseVariant(fields)) return;
  const int numIndiv=fields.size()-9;

  // The usual case: a genotype is listed for all individuals
  if(numIndiv==numIndividuals) 
    for(int i=0 ; i<numIndiv ; ++i) {
      const String &genotype=fields[i+9];
      Individual &indiv=individuals[i];
      indiv.genotypes.push_back(genotype);
    }

  // Odd case: genotypes are given only for males (dbSNP does this)
  else if(numIndiv==males.size() && males.size()>0) {
    int nextMale=0;
    for(int i=0 ; i<numIndiv ; ++i) {
      const String &genotype=fields[i+9];
      for( ; nextMale<numIndividuals &&
	     !males.isMember(individuals[nextMale].id) ; ++nextMale);
      if(nextMale>=numIndividuals) throw "too many fields in VCF file";
      Individual &indiv=individuals[nextMale];
      indiv.genotypes.push_back(genotype);
      ++nextMale;
    }
  }
  else throw String("Number of fields in VCF line does not match the number of individuals: ")+line;
}



bool Application::keep(const String &chr,int pos)
{
  if(!regions.isDefined(chr)) return false;
  const Vector<Region> &regs=regions[chr];
  for(Vector<Region>::const_iterator cur=regs.begin(), end=regs.end() ;
      cur!=end ; ++cur) {
    if((*cur).contains(pos)) return true;
  }
  return false;
}



void Application::output(File &out)
{
  // Output list of variants
  const int numVariants=variants.size();
  for(int i=0 ; i<numVariants ; ++i) {
    const Variant &v=variants[i];
    out.print(v.asString());
    out.print(i+1<numVariants ? "\t" : "\n");
  }
  if(numVariants==0) out.print("\n");

  // Output individual genotypes
  for(Vector<Individual>::const_iterator cur=individuals.begin(), 
	end=individuals.end() ; cur!=end ; ++cur) {
    const Individual &indiv=*cur;
    if(!keepIndiv.isMember(indiv.id)) continue;
    if(indiv.genotypes.size()==0) { out.print(indiv.id+"\n"); continue; }
    out.print(indiv.id+"\t");
    
    const int numVariants=indiv.genotypes.size();
    for(int i=0 ; i<numVariants ; ++i) {
      out.print(indiv.genotypes[i]);
      out.print(i+1<numVariants ? "\t" : "\n");
    }
  }
}



void Application::preprocess(File &infile)
{
  while(!infile.eof()) {
    String line=infile.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    if(fields.size()>0) {
      if(fields[0]=="#CHROM") parseChromLine(fields);
      else if(fields[0][0]!='#') parseVariant(fields);
    }
    delete &fields;
  }
}



bool Application::variableSite(const Vector<String> &fields)
{
  bool seenZero=false, seenOne=false;
  for(int i=9 ; i<fields.size() ; ++i) {
    const String &field=fields[i];
    if(field.length()==3) {
      if(field[0]=='0' || field[3]=='0') seenZero=true;
      if(field[0]!='0' || field[3]!='0') seenOne=true;
    }
    else if(field.length()==1) {
      if(field[0]=='0') seenZero=true;
      if(field[0]!='0') seenOne=true;
    }
    else throw String("Bad field in VCF: ")+field;
  }
  return seenZero && seenOne;
}



bool Application::parseVariant(const Vector<String> &fields,
			       String &chr,int &pos,String &ref,
			       Vector<String> &alts,String &id)
{
  if(variableOnly && !variableSite(fields)) return false;
  if(fields[6]!="PASS" && fields[6]!=".") return false;
  if(fields.size()<10) return false;
  if(fields[8]!="GT") {
    cerr<<"WARNING: GT expected in field 9 of VCF file"<<endl;
    return false;
  }
  if(fields[3].contains("<") || fields[4].contains("<")) return false;
  chr=fields[0];
  if(prependChr) chr=String("chr")+chr;
  pos=fields[1].asInt()-1; // VCF files are 1-based
  if(wantFilter && !keep(chr,pos)) return false;
  id=fields[2];
  if(id==".") id=chr+"@"+pos;
  ref=fields[3];
  const String &alt=fields[4];
  if(!dnaRegex.match(ref)) return false;
  if(SNPsOnly && (ref.size()!=1 || alt.size()!=1)) return false;
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



bool Application::parseVariant(const Vector<String> &fields)
{
  String chr, ref, id;
  Vector<String> alt;
  int pos;
  if(!parseVariant(fields,chr,pos,ref,alt,id)) return false;
  variants.push_back(Variant(chr,pos,ref,alt,id));
  return true;
}



void Application::loadIndivList(const String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw filename+" : cannot open file";
  while(!is.eof()) {
    String line;
    line.getline(is);
    Vector<String> fields;
    line.getFields(fields);
    if(fields.size()>=1) keepIndiv.insert(fields[0]);
  }
  numIndividuals=keepIndiv.size();
}



void Application::loadGender(const String &filename)
{
  knowMales=true;
  ifstream is(filename.c_str());
  int n=0;
  while(!is.eof()) {
    String line; line.getline(is);
    line.trimWhitespace();
    Vector<String> fields;
    line.getFields(fields);
    if(fields.size()==0) continue;
    if(fields.size()!=2) 
      throw String("Abort: ")+line+" : invalid gender specification";
    if(fields[1]=="male") males.insert(fields[0]);
    ++n;
  }
  if(numIndividuals==0) numIndividuals=n;
}
