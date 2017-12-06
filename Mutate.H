/****************************************************************
 Mutate.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Mutate_H
#define INCL_Mutate_H
#include "BOOM/Interval.H"
#include "BOOM/Array2D.H"
#include "BOOM/RouletteWheel.H"
#include "BOOM/TempFilename.H"
#include "BOOM/FastaWriter.H"
#include "ACEplus.H"
#include "ContentSensors.H"
#include "Model.H"
using namespace std;
using namespace BOOM;

struct SpliceSite {
  int pos;
  SignalType type;
  SpliceSite(SignalType t,int pos) : type(t), pos(pos) {}
};

struct SubstMatrix {
  Array2D<float> matrix;
  String substrate;
  Interval interval;
  Array1D<RouletteWheel> multinomials;
  SubstMatrix() : matrix(4,4), multinomials(4) {}
  virtual void prepareMultinomials();
  Symbol sample(Symbol from);
};

struct AllenMatrix : public SubstMatrix {
  AllenMatrix();
  void load(const String &filename);
  Symbol sample(int pos,const Sequence &chrom);
  virtual void prepareMultinomials();
};

class Mutate : public ACEplus {
public:
  Mutate();
  virtual ~Mutate() {}
  virtual int main(int argc,char *argv[]);
protected:
  int totalMutations;
  int brokenSites;
  int denovoDonors, denovoAcceptors;
  ContentSensors contentSensors;
  Model model;
  int MAX_GENES;
  int MUTATIONS_PER_GENE;
  String matrixFile;
  String gffTempFile, refTempFile, altTempFile, outTempFile, essexTempFile;
  float MIN_EXON_SCORE, MIN_INTRON_SCORE;
  FastaWriter fastaWriter;
  String bin;
  bool USE_SPLICE_MODEL;
  int MAX_DISTANCE;
  //Map<String,Vector<SubstMatrix> > matrices;
  AllenMatrix matrix;
  virtual void parseCommandLine(const CommandLine &);
  virtual void processConfig(const String &filename);
  virtual ContentSensor *loadContentSensor(const String &label,
					   ConfigFile &);
  void loadInputs(const String &configFile,const String &refGffFile,
		  const String &refFasta);
  void getSpliceSites(const Vector<GffTranscript*> &,Vector<SpliceSite> &into,
		      const String &seq);
  void simulate(const String &seqStr,const Vector<SpliceSite> &,
		const AllenMatrix &,Vector<GffTranscript*> &);
  //const SubstMatrix &getMatrix(const String &substrate,int begin,int end);
  //void loadSubstMatrices();
  bool checkBroken_old(const Sequence &,const String &seqStr,
		   const Vector<SpliceSite> &,const int mutationPos);
  bool checkBroken(const Sequence &,const String &seqStr,
		   const Vector<SpliceSite> &,const int mutationPos);
  bool checkDenovo(const Sequence &,const String &newSeqStr,
		   const String &oldSeqStr,const int mutationPos,
		   const Vector<SpliceSite> &,Vector<GffTranscript*> &);
  bool isConsensus(SignalType,const String &);
  bool hitsExistingSite(const int mutationPos,const Vector<SpliceSite> &);
  bool createsConsensus(const String &oldSeq,const String &newSeq,
			int pos,SignalType &);
  void incrementDenovo(SignalType);
  bool createsSite(const Sequence &,const String &seqStr,
		   int consensusPos,SignalType);
  float scoreDenovo(const Sequence &seq,const String &newSeqStr,
		    const String &oldSeqStr,
		    Vector<GffTranscript*> &transcripts);
  float scoreDenovo(const Sequence &seq,const String &newSeqStr,
		    Vector<GffTranscript*> &transcripts,
		    SignalType,int pos,ContentType &);
  float loadScore(GffTranscript &refTrans,
		  Vector<GffTranscript*> &transcripts);
  void getCoords(GffTranscript &transcript,SignalType,int signalPos,
		 int &begin,int &end,ContentType &);
  ContentType getContentType(GffTranscript &,int pos,Interval &);
};

#endif

