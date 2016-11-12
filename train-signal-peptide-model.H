/****************************************************************
 train-signal-peptide-model.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TRAIN_SIG_PEP_MDL_H
#define INCL_TRAIN_SIG_PEP_MDL_H
#include <string>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/FastaReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Array1D.H"
#include "BOOM/Strand.H"


class PeptideExon
{
  int begin, end;
  BOOM::String sequence;
  friend class ExonComparator;
public:
  PeptideExon(int b,int e);
  void loadSequence(const BOOM::String &substrate);
  BOOM::String &getSequence();
};



class ExonComparator : public BOOM::Comparator<PeptideExon>
{
  virtual bool equal(PeptideExon &a,PeptideExon &b) 
    {return a.begin==b.begin;}
  virtual bool greater(PeptideExon &a,PeptideExon &b)
    {return a.begin>b.begin;}
  virtual bool less(PeptideExon &a,PeptideExon &b)
    {return a.begin<b.begin;}
};



class SignalPeptide
{
  BOOM::Vector<PeptideExon> exons;
  BOOM::String substrate;
  Strand strand;
public:
  SignalPeptide(BOOM::String substrate,Strand);
  void addExon(int begin,int end);
  void sortExons();
  void loadSequence(const BOOM::String &substrate);
  BOOM::String getSequence();
};



struct Field
{
  BOOM::Map<char,float> aminoAcidFreqs;
  int fieldLength; // in acids
  int begin, end;  // in nucleotides
};



class Application
{
  BOOM::Map<BOOM::String,BOOM::Set<int> > transcriptsOnSubstrate;
  BOOM::Map<int,SignalPeptide*> peptides; // transcript ID -> SignalPeptide
  BOOM::Vector<BOOM::GffTranscript*> *transcripts;
  BOOM::Map<int,BOOM::GffTranscript*> transcriptsById;
  BOOM::Map<char,BOOM::Map<BOOM::String,float> > codonFreqs; // acid x codon->flt
  BOOM::Array1D<Field*> fields;
  int numFields;
  int numCodons;

  void initFields(const BOOM::String &);
  void loadSignalPeptideCoords(const BOOM::String &gffFile,
			       const BOOM::String &fastaFile);
  void loadTranscriptCoords(const BOOM::String &gffFile);
  void studyCodonUsage(const BOOM::String &filename);
  void updateCodonFreqs(const BOOM::String &transcript);
  void normalizeCodonFreqs();
  void studyAminoUsage();
  void studyAminoUsage(Field &,bool lastField);
  void generateModel(const BOOM::String &filename,
		     const BOOM::String &startCodonModelFile);
public:
  Application();
  int main(int argc,char *argv[]);
};


#endif
