/****************************************************************
 gene-LD.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GCF.H"
#include "BOOM/Map.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;

struct Individual {
  int hap[2][2]; // [which chrom][which gene]
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
  void recode(GCF &,int whichGene);
protected:
  Map<String,int> haploMap[2]; // for two genes
  Vector<String> haplotypes[2]; // for two genes
  Array1D<float> hapFreqs[2]; // for two genes
  Array1D<Individual> individuals;
  String genotypeToString(const Array1D<int> &genotype);
  void recode(int indivID,const GcfIndividual &,int whichChrom,int whichGene);
  float computeMI(int &numHaps1,int &numHaps2);
  void getHapFreqs(int whichGene);
  float entropy(int whichGene);
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
  CommandLine cmd(argc,argv,"q");
  if(cmd.numArgs()!=2)
    throw String("gene-LD <gene1.gcf> <gene2.gcf>\n\
  -q = quiet (no header)\n");
  const String file1=cmd.arg(0);
  const String file2=cmd.arg(1);
  const bool quiet=cmd.option('q');

  // Load files
  GCF gcf1(file1), gcf2(file2);
  const int numVar1=gcf1.numVariants(), numVar2=gcf2.numVariants();
  const int numIndiv=gcf1.numIndividuals();
  if(gcf2.numIndividuals()!=numIndiv) throw "files have different numbers of individuals";
  individuals.resize(numIndiv);

  // Re-code the haplotypes as integers (separately for each gene)
  recode(gcf1,0); recode(gcf2,1);

  // Compute mutual information
  float MI;
  int numHaps1, numHaps2;
  MI=computeMI(numHaps1,numHaps2);
  if(!isFinite(MI)) MI=0;

  // Generate output
  if(!quiet) cout<<"MI\t#variants1\t#variants2\t#haplo1\t#haplo2\t"<<endl;
  cout<<MI<<"\t"<<numVar1<<"\t"<<numVar2<<"\t"<<numHaps1<<"\t"<<numHaps2<<"\t"<<endl;

  return 0;
}



void Application::recode(GCF &gcf,int whichGene)
{
  const int numIndiv=gcf.numIndividuals();
  for(int i=0 ; i<numIndiv ; ++i) {
    const GcfIndividual &indiv=gcf.getIthIndividual(i);
    recode(i,indiv,0,whichGene);
    recode(i,indiv,1,whichGene);
  }
}



void Application::recode(int indivID,const GcfIndividual &indiv,int whichChrom,
			 int whichGene)
{
  String s=genotypeToString(indiv.chrom[whichChrom]);
  int hapID;
  if(!haploMap[whichGene].isDefined(s)) {
    hapID=haploMap[whichGene][s]=haplotypes[whichGene].size();
    haplotypes[whichGene].push_back(s);
  }
  else hapID=haploMap[whichGene][s];
  individuals[indivID].hap[whichChrom][whichGene]=hapID;
}



String Application::genotypeToString(const Array1D<int> &gt)
{
  const int n=gt.size();
  char str[n+1]; str[n]='\0';
  for(int i=0 ; i<n ; ++i) str[i]='0'+gt[i];
  return str;
}



float Application::computeMI(int &numHaps1,int &numHaps2)
{
  // Tabulate counts of invidividual haplotypes in each gene
  getHapFreqs(0); getHapFreqs(1);

  // Iterate over all combinations of haps from two genes
  numHaps1=haplotypes[0].size(); numHaps2=haplotypes[1].size();
  Map<String,int> comboCounts;
  const int numIndiv=individuals.size();
  for(int i=0 ; i<numIndiv ; ++i) {
    const Individual &indiv=individuals[i];
    const int chrom1_gene1=indiv.hap[0][0], chrom1_gene2=indiv.hap[0][1];
    const int chrom2_gene1=indiv.hap[1][0], chrom2_gene2=indiv.hap[1][1];
    const String combo1=String(chrom1_gene1)+" "+chrom1_gene2;
    const String combo2=String(chrom2_gene1)+" "+chrom2_gene2;
    if(!comboCounts.isDefined(combo1)) comboCounts[combo1]=0;
    if(!comboCounts.isDefined(combo2)) comboCounts[combo2]=0;
    ++comboCounts[combo1]; ++comboCounts[combo2];
  }
  Set<String> combinations;
  comboCounts.getKeys(combinations);
  float MI=0;
  for(Set<String>::iterator cur=combinations.begin(), end=combinations.end() ;
      cur!=end ; ++cur) {
    const String &combo=*cur;
    Vector<String> fields; combo.getFields(fields); if(fields.size()!=2) INTERNAL_ERROR;
    const int hap1=fields[0].asInt(), hap2=fields[1].asInt();
    const float pAB=comboCounts[combo]/float(numIndiv*2);
    const float pA=hapFreqs[0][hap1], pB=hapFreqs[1][hap2];
    MI+=pAB*log2(pAB/(pA*pB));
  }

  // Normalize by sqrt of product of entropies
  const float H1=entropy(0), H2=entropy(1);
  MI/=sqrt(H1*H2);
  return MI;
}



float Application::entropy(int whichGene)
{
  const int n=hapFreqs[whichGene].size();
  float H=0;
  for(int i=0 ; i<n ; ++i) {
    const float P=hapFreqs[whichGene][i];
    H-=P*log2(P);
  }
  return H;
}



void Application::getHapFreqs(int whichGene)
{
  const int numHaps=haplotypes[whichGene].size();
  hapFreqs[whichGene].resize(numHaps);
  hapFreqs[whichGene].setAllTo(0.0);
  const int numIndiv=individuals.size();
  for(int i=0 ; i<numIndiv ; ++i) {
    const Individual &indiv=individuals[i];
    ++hapFreqs[whichGene][indiv.hap[0][whichGene]]; // chrom 1
    ++hapFreqs[whichGene][indiv.hap[1][whichGene]]; // chrom 2
  }
  const int denom=numIndiv*2;
  for(int i=0 ; i<numHaps ; ++i) hapFreqs[whichGene][i]/=denom;
}

