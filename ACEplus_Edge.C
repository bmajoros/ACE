/****************************************************************
 ACEplus_Edge.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ACEplus_Edge.H"
#include "BOOM/Constants.H"
#include "LightVertex.H"
using namespace std;
using namespace BOOM;

ACEplus_Edge::ACEplus_Edge(const String &substrate,ContentType type,
		     LightVertex *from,LightVertex *to,int begin,int end,
		     Strand s,int ID)
  : LightEdge(substrate,type,from,to,begin,end,s,ID)
{
}



StructureChange &ACEplus_Edge::getChange()
{
  return change;
}



static String contentTypeString(ContentType type)
{
  switch(type) {
  case INTERGENIC: return "intergenic";
  case INTRON: return "intron";
  case EXON: return "exon";
  }
  INTERNAL_ERROR;
}



void ACEplus_Edge::printOn(ostream &os)
{
  int supportString=supported ? 1 : 0;
  bool useFrames=::isIntron(type) || ::isCoding(type);
  //bool needNewline=false;
  if(useFrames) {
    for(int i=0 ; i<3 ; ++i)
      if(isFinite(score[i])) {
	//if(needNewline) os<<endl;
	os<<substrate<<"\tedge\t"<<contentTypeString(type)<<"\t"<<begin+1
	  <<"\t"<<end<<"\t"<<score[i]<<"\t"<<strand<<"\t"<<i
	  <<"\tleft="<<left->getID()<<"; right="<<right->getID()
	  <<"; edgeID="<<ID<<"; sup="<<supportString<<";"
	  <<endl;
	//needNewline=true;
      }
  }
  else
    os<<substrate<<"\tedge\t"<<contentTypeString(type)<<"\t"<<begin+1
      <<"\t"<<end<<"\t"<<score[0]<<"\t"<<strand<<"\t.\tleft="
      <<left->getID()<<"; right="<<right->getID()<<"; edgeID="<<ID<<";"
      <<" sup="<<supportString<<";"
      <<" changes=\""<<change<<"\""
      <<endl;
}


