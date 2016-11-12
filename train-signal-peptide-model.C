/****************************************************************
 train-signal-peptide-model.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "train-signal-peptide-model.H"

Alphabet alphabet;

int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    BOOM::CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=5)
      throw BOOM::String(
"\ntrain-signal-peptide-model <*.gff> <*.fasta> <outfile>\n\
                           <start-codon.model> <field-lengths>\n\
   where <*.gff> contains items with the type \"signal_peptide\"\n\
         <window-len> = length in #codons\n\
         <start-codon.model> is a WMM\n\
         <field-lengths> is a comma-separate list of lengths (in acids)\n\
\n");
    BOOM::String gffFile=cmd.arg(0);
    BOOM::String fastaFile=cmd.arg(1);
    BOOM::String outfile=cmd.arg(2);
    BOOM::String startCodonModelFile=cmd.arg(3);
    BOOM::String lengths=cmd.arg(4);
    initFields(lengths);

    // Load signal peptide coordinates from GFF
    cerr<<"Loading signal peptide coordinates..."<<endl;
    loadSignalPeptideCoords(gffFile,fastaFile);

    // Load transcript coordinates from GFF
    cerr<<"Loading transcript coordinates..."<<endl;
    loadTranscriptCoords(gffFile);    
    
    // Study codon usage in transcripts
    cerr<<"Collecting codon usage statistics from transcripts..."<<endl;
    studyCodonUsage(fastaFile);
    
    // Study amino acid usage in signal peptides
    cerr<<"Collecting amino acid frequencies from signal peptides..."<<endl;
    studyAminoUsage();

    // Generate output
    cerr<<"Finishing model and writing output..."<<endl;
    generateModel(outfile,startCodonModelFile);

    return 0;
  }



void Application::loadSignalPeptideCoords(const BOOM::String &gffFile,
					  const BOOM::String &fastaFile)
{
  BOOM::GffReader reader(gffFile);
  BOOM::GffFeature *feature;
  while(feature=reader.nextFeature())
    {
      if(feature->getFeatureType()=="signal_peptide")
	{
	  Strand strand=feature->getStrand();
	  int transcriptId=feature->getTranscriptId();
	  if(!peptides.isDefined(transcriptId))
	    peptides[transcriptId]=
	      new SignalPeptide(feature->getSubstrate(),strand);
	  SignalPeptide *peptide=peptides[transcriptId];
	  int begin=feature->getBegin();
	  int end=feature->getEnd();
	  peptide->addExon(begin,end);
	  BOOM::String substrate=feature->getSubstrate();
	  transcriptsOnSubstrate[substrate].insert(transcriptId);
	}
      delete feature;
    }
}



void Application::loadTranscriptCoords(const BOOM::String &gffFile)
{
  BOOM::GffReader reader(gffFile);
  transcripts=reader.loadTranscripts();
  BOOM::Vector<BOOM::GffTranscript*> &trans=*transcripts;
  int n=transcripts->size();
  for(int i=0 ; i<n ; ++i)
    {
      BOOM::GffTranscript *transcript=trans[i];
      int transcriptId=transcript->getTranscriptId().asInt();
      transcriptsById[transcriptId]=transcript;
    }
}



void Application::studyCodonUsage(const BOOM::String &filename)
{
  BOOM::FastaReader reader(filename);
  BOOM::String defline, sequence;
  BOOM::Set<BOOM::String> substratesSeen;
  while(reader.nextSequence(defline,sequence))
    {
      BOOM::String substrateId, remainder;
      BOOM::FastaReader::parseDefline(defline,substrateId,remainder);
      if(substratesSeen.isMember(substrateId)) continue;
      substratesSeen.insert(substrateId);
      BOOM::Set<int> &transcriptIds=transcriptsOnSubstrate[substrateId];
      BOOM::Set<int>::iterator cur=transcriptIds.begin(), 
	end=transcriptIds.end();
      for(; cur!=end ; ++cur)
	{
	  int id=*cur;
	  BOOM::GffTranscript *transcript=transcriptsById[id];
	  transcript->loadSequence(sequence);
	  BOOM::String transSeq=transcript->getSequence();
	  updateCodonFreqs(transSeq);
	  SignalPeptide *peptide=peptides[id];
	  peptide->loadSequence(sequence);
	}
    }
  normalizeCodonFreqs();
}



void Application::updateCodonFreqs(const BOOM::String &transcript)
{
  const char *str=transcript.c_str();
  int len=transcript.length();
  const char *p=str;
  for(int i=0 ; i<len ; i+=3, p+=3)
    {
      BOOM::String codon(p,3);
      char acid=BOOM::ProteinTrans::mapCodon(codon.c_str());
      BOOM::Map<BOOM::String,float> &counts=codonFreqs[acid];
      if(!counts.isDefined(codon)) counts[codon]=1;
      else ++counts[codon];
    }
}



void Application::normalizeCodonFreqs()
{
  BOOM::Map<char,BOOM::Map<BOOM::String,float> >::iterator 
    cur=codonFreqs.begin(), end=codonFreqs.end();
  for(; cur!=end ; ++cur)
    {
      char acid=(*cur).first;
      BOOM::Map<BOOM::String,float> &theMap=(*cur).second;
      float total=0;
      BOOM::Map<BOOM::String,float>::iterator cur=theMap.begin(),
	end=theMap.end();
      for(; cur!=end ; ++cur)
	{
	  BOOM::String codon=(*cur).first;
	  float count=(*cur).second;
	  total+=count;
	}
      //cout<<acid<<" : ";
      cur=theMap.begin();
      for(; cur!=end ; ++cur)
	{
	  BOOM::String codon=(*cur).first;
	  float &count=(*cur).second;
	  count/=total;
	  //cout<<codon<<"="<<count<<" ";
	}
      //cout<<endl;
    }
}



void Application::studyAminoUsage()
{
  for(int i=0 ; i<numFields ; ++i)
    studyAminoUsage(*fields[i],i==numFields-1);
}



void Application::studyAminoUsage(Field &field,bool lastField)
{
  // Get field boundaries (in amino acid coordinates)
  int begin=field.begin/3, end=field.end/3;

  // Iterate through all peptides, counting the usage of the
  // individual amino acids in the respective fields
  int total=0;
  BOOM::Map<int,SignalPeptide*>::iterator pCur=peptides.begin(),
    pEnd=peptides.end();
  for(; pCur!=pEnd ; ++pCur)
    {
      SignalPeptide *peptide=(*pCur).second;
      const BOOM::String peptideSeq=peptide->getSequence();
      BOOM::String protein=BOOM::ProteinTrans::translate(peptideSeq);
      int len=protein.length();
      if(lastField) end=len;
      for(int i=begin ; i<end ; ++i)
	{
	  char acid=protein[i];
	  if(!field.aminoAcidFreqs.isDefined(acid))
	    field.aminoAcidFreqs[acid]=1;
	  else
	    field.aminoAcidFreqs[acid]+=1;
	  ++total;
	}
    }

  // Normalize frequencies into probabilities
  BOOM::Map<char,float>::iterator aCur=field.aminoAcidFreqs.begin(),
    aEnd=field.aminoAcidFreqs.end();
  for(; aCur!=aEnd ; ++aCur)
    {
      char acid=(*aCur).first;
      float &freq=(*aCur).second;
      freq/=total;
      //cout<<acid<<"="<<freq<<endl;
    }
}



void Application::generateModel(const BOOM::String &filename,
				const BOOM::String &startCodonModelFile)
{
  // Create output file
  ofstream os(filename.c_str());
  os<<"SignalPeptide"<<endl;

  // Copy the start codon model into the file
  ifstream is(startCodonModelFile.c_str());
  BOOM::String line;
  while(!is.eof())
    {
      line.getline(is);
      if(is.eof()) break;
      os<<line<<endl;
    }

  // Write out each field separately
  os<<numFields<<endl;
  for(int fieldNum=0 ; fieldNum<numFields ; ++fieldNum)
    {
      Field &field=*fields[fieldNum];

      // Count the number of codons (you just never know...)
      BOOM::Map<char,float>::iterator aCur=field.aminoAcidFreqs.begin(),
	aEnd=field.aminoAcidFreqs.end();
      numCodons=0;
      for(; aCur!=aEnd ; ++aCur)
	{
	  char acid=(*aCur).first;
	  float &acidP=(*aCur).second;
	  BOOM::Map<BOOM::String,float> codons=codonFreqs[acid];
	  numCodons+=codons.size();
	}

      os<<field.fieldLength<<endl;
      os<<numCodons<<endl;
      aCur=field.aminoAcidFreqs.begin();  aEnd=field.aminoAcidFreqs.end();
      for(; aCur!=aEnd ; ++aCur)
	{
	  char acid=(*aCur).first;
	  float &acidP=(*aCur).second;
	  BOOM::Map<BOOM::String,float> codons=codonFreqs[acid];
	  BOOM::Map<BOOM::String,float>::iterator cur=codons.begin(),
	    end=codons.end();
	  for(; cur!=end ; ++cur)
	    {
	      BOOM::String codon=(*cur).first;
	      float codonP=(*cur).second;
	      float logP=log(acidP*codonP);
	      os<<codon<<" "<<logP<<endl;
	    }
	}
    }
}



void Application::initFields(const BOOM::String &lengthString)
{
  BOOM::Vector<BOOM::String> &fieldLengths=*lengthString.getFields(",");
  numFields=fieldLengths.size();
  fields.resize(numFields);
  int begin=0;
  for(int i=0 ; i<numFields ; ++i)
    {
      int fieldLength=fieldLengths[i].asInt();
      Field *field=fields[i]=new Field;
      field->fieldLength=fieldLength;
      field->begin=begin;
      begin+=fieldLength*3;
      field->end=begin;
      //cout<<"field "<<i<<" : "<<fieldLength<<" "<<field->begin<<" "<<field->end<<endl;
    }
  delete &fieldLengths;
}



//===============================================================
//                    SignalPeptide methods
//===============================================================

SignalPeptide::SignalPeptide(BOOM::String substrate,Strand s) 
  : substrate(substrate), strand(s) 
{
  // ctor
}



void SignalPeptide::sortExons()
{
  ExonComparator cmp;
  BOOM::VectorSorter<PeptideExon> sorter(exons,cmp);
  switch(strand)
    {
    case FORWARD_STRAND:
      sorter.sortAscendInPlace();
      break;
    case REVERSE_STRAND:
      sorter.sortDescendInPlace();
      break;
    }
}



void SignalPeptide::addExon(int begin,int end)
{
  exons.push_back(PeptideExon(begin,end));
}



void SignalPeptide::loadSequence(const BOOM::String &substrate)
{
  int n=exons.size();
  for(int i=0 ; i<n ; ++i)
    exons[i].loadSequence(substrate);
}



BOOM::String SignalPeptide::getSequence()
{
  BOOM::String sequence;
  int n=exons.size();
  for(int i=0 ; i<n ; ++i)
    sequence+=exons[i].getSequence();
  if(sequence[0]=='M') 
    sequence=sequence.substring(1,sequence.length()-1);
  return sequence;
}



//===============================================================
//                    PeptideExon methods
//===============================================================
PeptideExon::PeptideExon(int b,int e) 
  : begin(b), end(e) 
{
  // ctor
}



void PeptideExon::loadSequence(const BOOM::String &substrate)
{
  sequence=substrate.substring(begin,end-begin);
}



BOOM::String &PeptideExon::getSequence()
{
  return sequence;
}


