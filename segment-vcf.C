/****************************************************************
 segment-vcf.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Pipe.H"
#include "BOOM/Regex.H"
#include "BOOM/VcfReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Interval.H"
#include "BOOM/BitSet.H"
using namespace std;
using namespace BOOM;

/*
class GffComparator : public Comparator<Interval> {
public:
  bool equal(Interval &a,Interval &b) { return a.getBegin()==b.getBegin(); }
  bool greater(Interval &a,Interval &b) { return a.getBegin()>b.getBegin(); }
  bool less(Interval &a,Interval &b) { return a.getBegin()<b.getBegin(); }
};
*/

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  bool haveGFF;
  BitSet bits;
  //Vector<Interval> features;
  String getChrom(const String &vcfFilename);
  void loadGFF(const String &filename,const String &chr);
  //void sortGFF();
  void emit(ostream &,const String &chr,int prevBoundary,int nextBoundary);
  void mask(int begin,int end);
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
  CommandLine cmd(argc,argv,"g:");
  if(cmd.numArgs()!=4)
    throw String("\n\
segment-vcf [options] <in.vcf> <binsize> <chrom-length> <out.bed>\n\
    -g file : also respect features given in GFF file\n\
\n\
    * VCF file must be for a single chromosome only\n\
    * GFF file can contain multiple chromosomes\n\
");
  //    * VCF file must be sorted by position\n\

  const String &infile=cmd.arg(0);
  const int binSize=cmd.arg(1).asInt();
  const int chromLen=cmd.arg(2).asInt();
  const String &outfile=cmd.arg(3);

  // Some initialization
  bits.setSize(chromLen);
  String chr=getChrom(infile);

  // Process optional GFF file
  haveGFF=cmd.option('g');
  if(haveGFF) {
    cout<<"loading GFF"<<endl; system("date");
    loadGFF(cmd.optParm('g'),chr);
    //sortGFF();
  }

  // Process VCF file
  /*
  ofstream os(outfile.c_str());
  VcfReader reader(infile);
  Variant prevVariant, nextVariant; Vector<Genotype> genotype;
  int featureIndex=0, prevBoundary=0, nextBoundary=0;
  reader.nextVariant(prevVariant,genotype);
  while(reader.nextVariant(nextVariant,genotype)) {
    nextBoundary=prevBoundary+binSize;
    if(nextBoundary>chromLen) { nextBoundary=chromLen; break; }
    if(nextVariant.getBegin()>nextBoundary) {
      if(nextBoundary<prevVariant.getEnd()) nextBoundary=prevVariant.getEnd();
      
    }
    prevVariant=nextVariant;
  }
  if(nextBoundary>prevBoundary)
    emit(os,chr,prevBoundary,nextBoundary);
  */
  cout<<"loading VCF"<<endl; system("date");
  VcfReader reader(infile);
  Variant v; Vector<Genotype> genotype;
  while(reader.nextVariant(v,genotype)) mask(v.getPos(),v.getEnd());

  // Segment
  cout<<"segmenting"<<endl; system("date");
  ofstream os(outfile.c_str());
  int prevPos=0;
  while(prevPos<chromLen) {
    int nextPos=prevPos+binSize;
    if(nextPos>=chromLen) nextPos=chromLen;
    while(nextPos<chromLen && bits.isMember(nextPos)) ++nextPos;
    emit(os,chr,prevPos,nextPos);
    prevPos=nextPos;
  }

  system("date");
  cout<<"[done]"<<endl;
  return 0;
}



String Application::getChrom(const String &vcfFilename)
{
  VcfReader reader(vcfFilename);
  Variant variant; Vector<Genotype> G;
  if(!reader.nextVariant(variant,G))
    throw "Error loading first variant from VCF file";
  return variant.getChr();
}



void Application::loadGFF(const String &filename,const String &chr)
{
  GffReader reader(filename);
  GffFeature *f;
  while(f=reader.nextFeature()) {
    if(f->getSubstrate()!=chr) continue;
    //features.push_back(Interval(f->getBegin(),f->getEnd()));
    mask(f->getBegin(),f->getEnd());
  }
}


/*
void Application::sortGFF()
{
  GffComparator cmp;
  VectorSorter<Interval> sorter(features,cmp);
  sorter.sortAscendInPlace();
}
*/


void Application::emit(ostream &os,const String &chr,int prevBoundary,
		       int nextBoundary)
{
  os<<chr<<"\t"<<prevBoundary<<"\t"<<nextBoundary<<endl;
}



void Application::mask(int begin,int end)
{
  if(begin>end) INTERNAL_ERROR;
  for(int i=begin ; i<end ; ++i) bits.addMember(i);
}





