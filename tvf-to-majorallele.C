/****************************************************************
 tvf-to-majorallele.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/Pipe.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Array1D.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/Regex.H"
#include "BOOM/Map.H"
#include "BOOM/TempFilename.H"
#include "Variant.H"
using namespace std;
using namespace BOOM;



inline int min(int a,int b) { return a<b ? a : b; }



// A Genotype represents the alleles of an individual at a single locus
struct Genotype {
  Array1D<int> alleles;
  Genotype(int ploidy) : alleles(ploidy) {}
  int ploidy() const { return alleles.size(); }
};




struct Region { // A chromosome or gene or other genomic interval
  String chr, id;
  int begin, end;
  String seq;
  char strand;
  Vector<Variant> variants;
  Region(const String &id,const String &chr,char strand,int begin,int end,
	 const String &seq) : id(id), chr(chr), begin(begin), end(end), 
			      strand(strand), seq(seq) {}
  bool contains(const Variant &v) const 
    { return v.chr==chr && v.refPos>=begin && v.refPos
	+v.alleles[0].getLength()<=end; }
  void loadSeq(const String &twoBitToFa,const String &genomeFile,
	       const String &tempFile);
  void clearSeq() { seq=""; }
  void printOn(ostream &os) const 
    { os<<chr<<":"<<begin<<"-"<<end<<":"<<strand; }
};
ostream &operator<<(ostream &os,const Region &r) { r.printOn(os); return os; }



bool operator<(const Variant &v,const Region &r)
{ return v.refPos<r.begin; }
bool operator>(const Variant &v,const Region &r)
{ return v.refPos>r.end; }



// This assumes no overlaps, so we need not compare ends
struct RegionComp : Comparator<Region*> {
  bool equal(Region *&a,Region *&b)
  { return a->chr==b->chr && a->begin==b->begin; }
  bool greater(Region *&a,Region *&b) 
  { return a->chr>b->chr || a->chr==b->chr && a->begin>b->begin; }
  bool less(Region *&a,Region *&b)    
  { return a->chr<b->chr || a->chr==b->chr && a->begin<b->begin; }
};



class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  bool SANITY_CHECKS, DRY_RUN;
  int PLOIDY; // default=2, can override on command line with -p
  String twoBitToFa;
  String tempfile;
  Regex gzRegex;
  Map<String,Vector<Region*> > regionsByChr;
  Vector<Region*> regions;
  Vector<Variant> variants;
  FastaWriter writer;
  bool wantRef;
  bool nonhuman;
  String wantIndiv, wantChr, genomeFile;
  Set<String> males;
  bool knowMales; // whether the list of males was given
  void loadMales(const String &);
  void convert(File &tvf,ostream *,const String genomeFile);
  void parseHeader(const String &line);
  void loadRegions(const String &regionsFilename,const String &genomeFilename,
		   const String &tempFile);
  void emit(const String &individualID,const Vector<Genotype> &loci,ostream *);
  void updateCigar(int refLen,int altLen,int localPos,int &matchBegin,
		   String &cigar,int delta);
  bool skipGender(const Region &region,bool male,int ploid);
  const Variant *disambiguateOverlaps(int &v,const int numVariants,int ploid,
				      const Vector<Variant> &variants,
				      const Vector<Genotype> &loci,
				      int delta,const String &altGenome,
				      int regionBegin,const String &indiv,
				      const String &regionID,
				      bool &refMismatch);
  void transitiveClosure(const Vector<Variant> &variants,int &v,
			 const int numVariants,int ploid,
			 const Vector<Genotype> &loci,
			 Vector<const Variant*> &closure);
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
  : twoBitToFa("twoBitToFa"), gzRegex("gz$"), tempfile(TempFilename::get()),
    PLOIDY(2), SANITY_CHECKS(true), DRY_RUN(false)
{
  // ctor
  randomize();
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"dhrt:i:c:y:p:s");
  if(cmd.numArgs()!=5)
    throw String("\ntvf-to-fasta [options] <in.tvf> <genome.2bit> <regions.bed> <out.fasta> <sample-info.txt>\n\
     sample-info.txt: tab-separated file:\n\
          indiv gender subpopulation superpopulation\n\
     -d : dry run - no output, just report errors\n\
     -t path : path to twoBitToFa\n\
     -r : emit reference sequence also\n\
     -i ID : only this sample (individual)\n\
     -c chr : only regions on this chromosome\n\
     -y gender.txt : one individual per line with \"male\" or \"female\"\n\
     -h : nonhuman species (don't treat chrX/chrY by gender)\n\
     -p ploidy : default is 2\n\
     -s : no sanity checks (improve speed)\n\
\n\
     NOTE: Run 'which twoBitToFa' to ensure it's in your path\n\
     NOTE: regions.bed is a BED6 file: chr begin end name score strand\n\
     NOTE: For human subjects, chrom name must begin with \"chr\"\n\
");
  const String &tvfFilename=cmd.arg(0);
  genomeFile=cmd.arg(1);
  const String &regionsFilename=cmd.arg(2);
  const String &fastaFilename=cmd.arg(3);
  const String &sampleInfoFile=cmd.arg(4);
  if(cmd.option('t')) twoBitToFa=cmd.optParm('t');
  wantRef=cmd.option('r');
  if(cmd.option('i')) wantIndiv=cmd.optParm('i');
  if(cmd.option('c')) wantChr=cmd.optParm('c');
  if(cmd.option('y')) loadMales(cmd.optParm('y'));
  nonhuman=cmd.option('h');
  if(cmd.option('p')) PLOIDY=cmd.optParm('p').asInt();
  if(cmd.option('s')) SANITY_CHECKS=false;
  if(cmd.option('d')) DRY_RUN=true;

  // Load regions
  loadRegions(regionsFilename,genomeFile,fastaFilename);

  // Process TVF file
  File *tvf=gzRegex.search(tvfFilename) ? new GunzipPipe(tvfFilename)
    : new File(tvfFilename);
  ofstream *os=DRY_RUN ? NULL : new ofstream(fastaFilename.c_str());
  convert(*tvf,os,genomeFile);
  delete tvf;
  delete os;

  return 0;
}



void Application::convert(File &tvf,ostream *os,const String genomeFile)
{
  // Parse header
  String line=tvf.getline();
  line.trimWhitespace();
  parseHeader(line);
  const int numVariants=variants.size();

  // Emit reference
  if(wantRef) {
    for(Vector<Region*>::const_iterator cur=regions.begin(), 
	  end=regions.end() ; cur!=end ; ++cur) {
      const Region &region=**cur;
      String seq=region.seq;
      //if(region.strand=='-') seq=ProteinTrans::reverseComplement(seq);
      const int L=seq.length();
      const String cigar=String("")+L+"M";
      for(int j=0 ; j<PLOIDY ; ++j) {
	int alleleNum=j+1;
	String def=String(">reference_")+j+" /individual=reference"+
	  " /allele="+alleleNum+" /locus="+region.id+" /coord="+region.chr+":"
	  +region.begin+"-"+region.end+":"+region.strand+" /cigar="+cigar
	  +" /variants=";
	if(os) writer.addToFasta(def,seq,*os);
      }
    }
  }

  // Process each individual
  while(!tvf.eof()) {
    line=tvf.getline();
    line.trimWhitespace();
    if(line.isEmpty()) continue;
    Vector<String> &fields=*line.getFields();
    if(fields.size()!=numVariants+1) {
      cout<<fields.size()<<"\t"<<numVariants<<"\t"<<endl;
      if(fields.size()>0) cout<<fields[0]<<endl;
      throw String("Abort: Wrong number of fields in TVF line: ")+line;
    }
    String id=fields.front();
    if(!wantIndiv.isEmpty() && id!=wantIndiv) continue;
    fields.erase(fields.begin());
    Vector<Genotype> loci;
    for(Vector<String>::const_iterator cur=fields.begin(), end=fields.end() ;
	cur!=end ; ++cur) {
      const String &field=*cur;
      if(field.getLength()==1) {
	Genotype gt(1);
	if(field[0]=='.') gt.alleles[0]=0;
	else {
	  gt.alleles[0]=field[0]-'0';
	  if(gt.alleles[0]<0 || gt.alleles[0]>9) 
	    throw String("Abort: ")+fields[0]+" : unknown allele indicator";
	}
	loci.push_back(gt);
      }
      else if(field.getLength()==3) { // ### ASSUMES DIPLOID
	if(field[1]!='|') throw "Abort: VCF file is not phased";
	Genotype gt(2);
	if(field[0]=='.') gt.alleles[0]=0;
	else { 
	  gt.alleles[0]=field[0]-'0'; 
	  if(gt.alleles[0]<0 || gt.alleles[0]>9) 
	    throw String("Abort: ")+fields[0]+" : unknown allele indicator";
	}
	if(field[2]=='.') gt.alleles[1]=0;
	else {
	  gt.alleles[1]=field[2]-'0';
	  if(gt.alleles[1]<0 || gt.alleles[1]>9) 
	    throw String("Abort: ")+fields[2]+" : unknown allele indicator";
	}
	loci.push_back(gt);
      }
      else throw String("Abort: Cannot parse genotype: ")+field;
    }
    delete &fields;
    emit(id,loci,os);
    if(!wantIndiv.isEmpty()) break;
  }
}



void Application::parseHeader(const String &line)
{
  Vector<String> fields; line.getFields(fields);
  Map<String,int> prevPos;
  int i=0;
  for(Vector<String>::const_iterator cur=fields.begin(), end=fields.end()
	; cur!=end ; ++cur, ++i) {
    const String &field=*cur;
    Vector<String> fields; field.getFields(fields,":");
    if(fields.size()<5) throw String("Abort: cannot parse TVF header: ")+field;
    const String &id=fields[0];
    const String &chr=fields[1];
    //if(!wantChr.isEmpty() && chr!=wantChr) continue;
    const int pos=fields[2];
    if(!prevPos.isDefined(chr)) prevPos[chr]=0;
    //if(pos==prevPos[chr]) continue;
    if(pos<prevPos[chr]) 
      throw String("Abort: ")+String(pos)+"<="+prevPos[chr]+
	": input file is not sorted: use vcf-sort and re-convert to tvf";
    prevPos[chr]=pos;
    Variant variant(id,chr,pos,0,i);
    for(int j=3 ; j<fields.size() ; ++j) variant.addAllele(fields[j]);
    variant.trim();
    variants.push_back(variant);
  }

  // Now add to regions
  for(Vector<Variant>::iterator vcur=variants.begin(), 
	vend=variants.end() ; vcur!=vend ; ++vcur) {
    const Variant &v=*vcur;
    const int pos=v.refPos;
    if(!regionsByChr.isDefined(v.chr)) continue;
    Vector<Region*> &rs=regionsByChr[v.chr];
    const int N=rs.size();
    for(int b=0, e=N ; b<e ; ) {
      const int mid=(b+e)/2;
      Region &midRegion=*rs[mid];
      if(pos<midRegion.begin) e=mid;
      else if(pos>=midRegion.end) b=mid+1;
      else { midRegion.variants.push_back(v); break; }
    }
  }
}



void Region::loadSeq(const String &twoBitToFa,const String &genomeFile,
		     const String &tempFile)
{
  // Invoke twoBitToFa to extract sequence from chrom file
  String cmd=twoBitToFa+" -seq="+chr+" -start="+begin+" -end="+end+
    +" "+genomeFile+" "+tempFile;
  system(cmd.c_str());
  String def;
  FastaReader::load(tempFile,def,seq);
}



void Application::loadRegions(const String &regionsFilename,const String &
			      genomeFilename,const String &tempFile)
{
  File reg(regionsFilename);
  while(!reg.eof()) {
    String line=reg.getline();
    line.trimWhitespace();
    if(line.isEmpty()) continue;
    Vector<String> fields; line.getFields(fields);
    if(fields.size()<6) throw "Abort: regions file requires 6 fields: chr begin end name score strand";
    const String chr=fields[0];
    if(!wantChr.isEmpty() && chr!=wantChr) continue;
    const int begin=fields[1].asInt(), end=fields[2].asInt();
    const String id=fields[3];
    char strand=fields[5][0];
    
    // Invoke twoBitToFa to extract sequence from chrom file
    String seq;
    if(wantIndiv.isEmpty()) {
      String cmd=twoBitToFa+" -seq="+chr+" -start="+begin+" -end="+end+
	+" "+genomeFilename+" "+tempFile;
      system(cmd.c_str());
      String def;
      FastaReader::load(tempFile,def,seq);
      //if(seq.contains(",")) throw String("Abort: ")+tempFile+"contains a comma: "+cmd;
    }

    Region *r=new Region(id,chr,strand,begin,end,seq);
    regions.push_back(r);
    regionsByChr[chr].push_back(r);
  }
  unlink(tempFile.c_str());

  // Sort regions so we can do fast searches later
  RegionComp cmp;
  Set<String> keys; regionsByChr.getKeys(keys);
  for(Set<String>::const_iterator cur=keys.begin(),
	end=keys.end() ; cur!=end ; ++cur) {
    VectorSorter<Region*> sorter(regionsByChr[*cur],cmp);
    sorter.sortAscendInPlace();
  }
}



void Application::loadMales(const String &filename)
{
  knowMales=true;
  ifstream is(filename.c_str());
  while(!is.eof()) {
    String line; line.getline(is);
    line.trimWhitespace();
    Vector<String> fields;
    line.getFields(fields);
    if(fields.size()==0) continue;
    if(fields.size()!=2) 
      throw String("Abort: ")+line+" : invalid gender specification";
    if(fields[1]=="male") males.insert(fields[0]);
  }
}



/*
  This function fins the transitive closure of the overlap relation, but
  only among those variants for which this individual has an alternate allele.
  It also advances the v index, so that upon return v will either be a valid
  index of the next alternate-allele variant to consider, or will be at
  the end of the array of variants.  Either way, v is always incremented at
  least once.  We assume variants are sorted by increasing position.
 */
void Application::transitiveClosure(const Vector<Variant> &variants,int &v,
				    const int numVariants,int ploid,
				    const Vector<Genotype> &loci,
				    Vector<const Variant*> &closure)
{
  int end=0;
  for( ; v<numVariants ; ++v) {
    const Variant &variant=variants[v];
    const int state=loci[variant.i].alleles[ploid];
    if(state==0) continue; // same as ref, nothing to do
    if(end>0 && variant.refPos>=end) break;
    const int refAlleleEnd=variant.refEnd();
    if(refAlleleEnd>end) end=refAlleleEnd;
    closure.push_back(&variant);
  }
}



/*
  This function disambiguates the case of one or more variants overlapping
  in the reference.  If all such overlapping variants are consistent, the 
  longest variant (by reference length) is returned and applied, which should
  be consistent with all enclosed variants.  If inconsistent variants are 
  detected, an error is reported.  Multiple insertions at the same location, 
  or an insertion within a deletion, or two SNPs at the same position are not
  allowed, and are reported as an error.  Partially overlapping deletions (in
  which neither is completely contained within the other) are rejected.
  In all cases, variants for which this individual has the ref allele are 
  ignored.
 */
const Variant *Application::disambiguateOverlaps(int &v,const int numVariants,
					 int ploid,
					 const Vector<Variant> &variants,
					 const Vector<Genotype> &loci,
					 int delta,const String &altGenome,
					 int regionBegin,const String &indiv,
					 const String &regionID,
					 bool &refMismatch)
{
  refMismatch=false;
  const int altGenomeLen=altGenome.length();

  // First, take the transitive closure of overlapping variants
  Vector<const Variant*> closure;
  transitiveClosure(variants,v,numVariants,ploid,loci,closure);
  if(closure.isEmpty()) return NULL;

  // Find the longest variant (by ref length) in the closure (i.e., the
  // longest deletion)
  int longestLength=-1;
  Variant *longestVariant=NULL;
  for(Vector<const Variant*>::iterator cur=closure.begin(), end=closure.end() ;
      cur!=end ; ++cur) {
    const Variant &variant=**cur;
    const int L=variant.refLen();
    if(L>longestLength) {
      longestLength=L;
      longestVariant=&variant; }
  }
  if(longestLength<0) INTERNAL_ERROR;

  // If all variants have length <=1 in the reference, take the longest
  // by alt length (i.e., the longest insertion)
  if(longestLength<=1) {
    longestVariant=NULL;
    longestLength=-1;
    for(Vector<const Variant*>::iterator cur=closure.begin(), 
	  end=closure.end() ; cur!=end ; ++cur) {
      const Variant &variant=**cur;
      const int state=loci[variant.i].alleles[ploid];
      const int L=variant.alleles[state].length();
      if(L>longestLength) {
	longestLength=L;
	longestVariant=&variant; }
    }
  }

  if(SANITY_CHECKS) {
    // Verify that all ref alleles match genomic sequence
    for(Vector<const Variant*>::iterator cur=closure.begin(), 
	  end=closure.end() ; cur!=end ; ++cur) {
      const Variant &variant=**cur;
      String ref=variant.alleles[0];
      int refLen=ref.length();
      const int altPos=variant.refPos-regionBegin-delta;
      if(altPos+refLen>altGenomeLen)
	{ refLen=altGenomeLen-altPos; ref=ref.substring(0,refLen); }
      String genomic=altGenome.substring(altPos,refLen);
      if(ref!=genomic) {
	cerr<<"VCF_ERROR\tSEQUENCE_MISMATCH\t"<<indiv<<"\t"<<regionID
	    <<"\t"<<variant.id<<":"<<variant.chr<<":"<<variant.refPos
	    <<"\t"<<ref<<"!="<<genomic<<endl;
	refMismatch=true;
	return NULL;
      }
    }

    // Verify that all others are strictly contained within the longest
    for(Vector<const Variant*>::iterator cur=closure.begin(), 
	  end=closure.end() ; cur!=end ; ++cur) {
      const Variant &variant=**cur;
      const int state=loci[variant.i].alleles[ploid];
      for(Vector<const Variant*>::iterator cur=closure.begin(), 
	    end=closure.end() ; cur!=end ; ++cur) {
	const Variant &other=**cur;
	if(&other==&variant) continue;
	if(variant.overlaps(other)) {
	  const int otherState=loci[other.i].alleles[ploid];
	  if(!other.covers(variant) && !variant.covers(other)) {
	    cerr<<"VCF_ERROR\tPARTIALLY_OVERLAPPING_VARIANTS\t"
		<<indiv<<"\t"<<regionID<<"\t";
	    variant.printAllele(state,cerr); cerr<<"\t";
	    other.printAllele(otherState,cerr); cerr<<endl;
	    return NULL; 
	  }
	}
      }
    }
    
    // Verify that nested variants are compatible
    for(Vector<const Variant*>::iterator cur=closure.begin(), 
	  end=closure.end() ; cur!=end ; ++cur) {
      const Variant &variant=**cur;
      const int state=loci[variant.i].alleles[ploid];
      for(Vector<const Variant*>::iterator cur=closure.begin(), 
	    end=closure.end() ; cur!=end ; ++cur) {
	const Variant &other=**cur;
	if(&other==&variant) continue;
	if(other.overlaps(variant)) {
	  const int otherState=loci[other.i].alleles[ploid];
	  if(variant.identical(state,other,otherState)) continue;
	  if(other.covers(variant) &&
	     other.alleles[otherState].contains(variant.alleles[state]) ||
	     variant.covers(other) &&
	     variant.alleles[state].contains(other.alleles[otherState])) {
	    cerr<<"VCF_WARNING\tNESTED_VARIANTS\t"<<indiv<<"\t"
		<<regionID<<"\t";
	    variant.printAllele(state,cerr); cerr<<"\t";
	    other.printAllele(otherState,cerr); cerr<<endl;
	    continue; // They are compatible
	  }
	  cerr<<"VCF_ERROR\tINCOMPATIBLE_NESTED_VARIANTS\t"<<indiv<<"\t"
	      <<regionID<<"\t";
	  variant.printAllele(state,cerr); cerr<<"\t";
	  other.printAllele(otherState,cerr); cerr<<endl;
	  return NULL; // Incompatible
	}
      }
    }
    
    // Verify that there are no insertions at the same position as another 
    // insertion (unless they're consistent)
    for(Vector<const Variant*>::iterator cur=closure.begin(), 
	  end=closure.end() ; cur!=end ; ++cur) {
      const Variant &variant=**cur;
      const int state=loci[variant.i].alleles[ploid];
      if(!variant.insertion(state)) continue;
      for(Vector<const Variant*>::iterator cur=closure.begin(), 
	    end=closure.end() ; cur!=end ; ++cur) {
	const Variant &other=**cur;
	if(&other==&variant) continue;
	if(variant.refPos==other.refPos) {
	  const int otherState=loci[other.i].alleles[ploid];
	  // Insertion overlapping another insertion:
	  if(other.insertion(otherState)) {
	    if(variant.alleles[state].contains(other.alleles[otherState]) ||
	       other.alleles[otherState].contains(variant.alleles[state])) {
	      cerr<<"VCF_WARNING\tNESTED_INSERTIONS\t"<<indiv<<"\t"
		  <<regionID<<"\t";
	      variant.printAllele(state,cerr); cerr<<"\t";
	      other.printAllele(otherState,cerr); cerr<<endl;	      
	      continue; // They are compatible
	    }
	    else {
	      cerr<<"VCF_ERROR\tAMBIGUOUS_INSERTIONS\t"<<indiv<<"\t"
		  <<regionID<<"\t";
	      variant.printAllele(state,cerr); cerr<<"\t";
	      other.printAllele(otherState,cerr); cerr<<endl;
	      return NULL; // Incompatible
	    }
	  }
	}
      }
    }
  }

  // Return the longest variant
  return longestVariant;
}



bool Application::skipGender(const Region &region,bool male,int ploid)
{
  if(nonhuman) return false;
  bool chrX=region.chr=="chrX", chrY=region.chr=="chrY";
  int ploidy=PLOIDY;
  if(male) { if(chrX || chrY) ploidy=1; }
  else { // female
    if(chrX) ploidy=2;
    else if(chrY) ploidy=0; }
  return ploid>=ploidy;
}



void Application::updateCigar(int refLen,int altLen,int localPos,
			      int &matchBegin,String &cigar,int delta)
{
  // Precondition: delta is nonzero
  int minAllele=min(refLen,altLen);
  const int matchEnd=localPos+minAllele;
  const int matchLen=matchEnd-matchBegin;
  if(matchLen>0) cigar+=String("")+matchLen+"M";
  const char indel=refLen>altLen ? 'D' : 'I';
  const int indelLen=delta>0 ? delta : -delta;
  cigar+=String("")+indelLen+indel;
  matchBegin=localPos+refLen;
}



void Application::emit(const String &individualID,const Vector<Genotype> &loci,
		       ostream *os)
{
  bool male=knowMales ? males.isMember(individualID) : true;

  // Handle each genome in this individual separately
  for(int ploid=0 ; ploid<PLOIDY ; ++ploid) {

    // Iterate over regions
    for(Vector<Region*>::const_iterator cur=regions.begin(), 
	  end=regions.end() ; cur!=end ; ++cur) {
      const Region &region=**cur;
      if(skipGender(region,male,ploid)) continue;
      int deltas=0;
      String cigar;
      if(!wantIndiv.isEmpty()) region.loadSeq(twoBitToFa,genomeFile,tempfile);
      unlink(tempfile.c_str());
      String seq=region.seq;
      int matchBegin=0;
      String region_hap=region.id+"_"+ploid;

      // Iterate over variants
      const int numVariants=region.variants.size();
      int variantsApplied=0, indelVariantsApplied=0, mismatches=0;
      String deflineVariants;
      for(int v=0 ; v<numVariants ; ) {
	bool refMismatch;
	const Variant *variant=
	  disambiguateOverlaps(v,numVariants,ploid,region.variants,loci,deltas,
			       seq,region.begin,individualID,region_hap,
			       refMismatch);
	if(refMismatch) ++mismatches;
	if(!variant) continue;

	// Prepare to do the substitution
	const int localPos=variant->refPos-region.begin;
	const int state=loci[variant->i].alleles[ploid];
	if(SANITY_CHECKS) if(!state) INTERNAL_ERROR;
	const String &refAllele=variant->alleles[0];
	const String &altAllele=variant->alleles[state];
	int refLen=refAllele.getLength(), altLen=altAllele.getLength();
	if(localPos+refLen>region.seq.getLength()) 
	  refLen=region.seq.getLength()-localPos;

	// Do the substitution in the alt genome
	seq.replaceSubstring(localPos-deltas,refLen,altAllele);
	++variantsApplied;

	// Add to the defline
	if(!deflineVariants.empty()) deflineVariants+=",";
	deflineVariants+=variant->id+":"+variant->chr+":"+localPos
	  +":"+(localPos-deltas)+":"+refAllele+":"+altAllele;

	// Update the delta (difference in coordinates btwn ref/alt)
	const int delta=refLen-altLen;
	deltas+=delta;
	if(delta>0) ++indelVariantsApplied;

	// Update CIGAR string
	if(delta!=0)
	  updateCigar(refLen,altLen,localPos,matchBegin,cigar,delta);
      } // end of foreach variant

      // Final update to CIGAR string
      const int matchLen=region.seq.length()-matchBegin;
      if(matchLen>0) cigar+=String("")+matchLen+"M";

      // Write into FASTA file
      String def=String(">")+individualID+"_"+ploid+" /individual="+
	individualID+" /allele="+(ploid+1)+" /locus="+region.id+" /coord="+
	region.chr+":"+region.begin+"-"+region.end+":"+
	region.strand+" /cigar="+cigar+" /variants="+deflineVariants;
      if(os) writer.addToFasta(def,seq,*os);
      if(!wantIndiv.isEmpty()) region.clearSeq(); // save memory

      // Report stats
      if(variantsApplied>0 || indelVariantsApplied>0 || mismatches>0)
	cout<<individualID<<"\thap"<<ploid<<"\t"<<region.id<<"\t"
	    <<variantsApplied<<" applied\t"
	    <<indelVariantsApplied<<" indels applied\t"<<mismatches
	    <<" mismatches"<<endl;
    } // end foreach region
  } // end for ploidy
}



/************************************************************************
  This function trims the ref and alt alleles so that they have no
  common prefix nor suffix.  This may result in an allele having zero
  length, which corresponds to a "pure" insertion or deletion.  Trimming
  is necessary to ensure that we can detect when variants that appear
  to overlap are actually consistent.
 */
void Variant::trim()
{
  //bool debug=alleles[0].length()>1 || alleles[1].length()>1;
  //if(debug) cout<<"BEFORE: "<<*this<<endl;

  const int numAlleles=alleles.size();
  
  // First, identify common sequence on the right
  const String &ref=alleles[0];
  Array1D<int> end(numAlleles);
  for(int i=0 ; i<numAlleles ; ++i) end[i]=alleles[i].length()-1;
  while(1) {
    bool stop=false;
    for(int i=0 ; i<numAlleles ; ++i) if(end[i]<0) { stop=true; break; }
    if(stop) break;
    const char r=ref[end[0]];
    bool match=true;
    for(int i=1 ; i<numAlleles ; ++i)
      if(alleles[i][end[i]]!=r) { match=false; break; }
    if(!match) break;
    for(int i=0 ; i<numAlleles ; ++i) --end[i];
  }
  // Now everything after end[i] is a common suffix
  //for(int i=0 ; i<numAlleles ; ++i) cout<<"end["<<i<<"]="<<end[i]<<endl;

  // Next, identify length of common sequence on the left
  int begin=0; // first non-matching position
  //const int refLen=ref.length();
  //while(begin<refLen) {
  while(1) {
    bool stop=false;
    for(int i=0 ; i<numAlleles ; ++i) if(begin>end[i]) { stop=true; break; }
    if(stop) break;
    const char r=ref[begin];
    bool match=true;
    for(int i=1 ; i<numAlleles ; ++i) {
      const String &alt=alleles[i];
      if(alt[begin]!=r) { match=false; break; }
    }
    if(!match) break;
    ++begin;
  }
  //cout<<"begin="<<begin<<endl;
  // begin now points to the first column in which either there is a 
  // mismatch or we're past the end of one allele (ref or alt)

  // Trim the alleles
  for(int i=0 ; i<numAlleles ; ++i) {
    int len=end[i]+1-begin;
    if(len<0) len=0;
    //cout<<"trimming "<<alleles[i];
    alleles[i]=alleles[i].substring(begin,len);
    //cout<<" to "<<alleles[i]<<endl;
  }
  //cout<<"================"<<endl;

  //for(int i=0 ; i<numAlleles ; ++i) cout<<"allele["<<i<<"]="<<alleles[i]<<endl;

  // Finally, adjust the position of this variant
  refPos+=begin;

  //if(debug) cout<<"AFTER: "<<*this<<endl<<endl;
}


