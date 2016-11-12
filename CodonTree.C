/****************************************************************
 CodonTree.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "CodonTree.H"
#include "BOOM/Constants.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/ProteinTrans.H"
#include <iostream>


CodonTree::CodonTree(istream &is)
  : array(125)
{
  load(is);
}



CodonTree::CodonTree(const CodonTree &other)
  : array(125)
{
  for(int i=0 ; i<125 ; ++i)
    array[i]=other.array[i];
}



void CodonTree::load(istream &is)
{
  for(int i=0 ; i<125 ; ++i)
    array[i]=NEGATIVE_INFINITY;
  int n;
  is>>n;
  char codon[4];
  codon[3]='\0';
  double logP;
  for(int i=0 ; i<n ; ++i)
    {
      is>>codon[0]>>codon[1]>>codon[2]>>logP;
      int index=codonToIndex(codon);
      array[index]=logP;
    }
}



int CodonTree::codonToIndex(const char *codon)
{
  Alphabet &alphabet=DnaAlphabet::global();
  int index=
    alphabet.lookup(codon[0])*25+
    alphabet.lookup(codon[1])*5+
    alphabet.lookup(codon[2]);
  return index;
}



float CodonTree::scoreCodon(const char *str,int pos)
{
  const char *codonStart=str+pos;
  int index=codonToIndex(codonStart);
  return array[index];
}



float CodonTree::scoreCodon(const char *str)
{
  return array[codonToIndex(str)];
}



void CodonTree::revComplementSelf()
{
  BOOM::Array1D<double> newArray(125);
  for(int index=0 ; index<125 ; ++index)
    {
      BOOM::String codon=decode(index);
      BOOM::String anticodon=BOOM::ProteinTrans::reverseComplement(codon);
      int newIndex=codonToIndex(anticodon.c_str());
      newArray[newIndex]=array[index];
    }
  for(int index=0 ; index<125 ; ++index)
    array[index]=newArray[index];
}



BOOM::String CodonTree::decode(int base4)
{
  Alphabet &alphabet=DnaAlphabet::global();
  BOOM::String codon="   ";
  int divisor=25;
  for(int i=0 ; i<3 ; ++i)
    {
      int x=base4/divisor;
      codon[i]=alphabet.lookup(x);
      base4-=x*divisor;
      divisor/=5;
    }
  return codon;
}



CodonTree *CodonTree::clone() const 
{
  return new CodonTree(*this);
}


