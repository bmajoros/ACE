/****************************************************************
 ProjectionChecker.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/ProteinTrans.H"
#include "BOOM/CodonIterator.H"
#include "BOOM/Sequence.H"
#include "BOOM/Essex.H"
#include "Labeling.H"
#include "TranscriptSignals.H"
#include "SignalSensors.H"
using namespace std;
using namespace BOOM;



class ProjectionChecker {
public:
  ProjectionChecker(const GffTranscript &refTrans,
		    const GffTranscript &altTrans,
		    const String &refSubstrate,
		    const Sequence &refSeq,
		    const String &altSubstrate,
		    const Sequence &altSeq,
		    const Labeling &labeling,
		    const SignalSensors &);
  bool checkSpliceSites(bool quiet);
  TranscriptSignals *findBrokenSpliceSites(); // NULL if bad mapping
  TranscriptSignals *simulateBrokenSpliceSites();
  bool checkSpliceSiteStrengths();
  bool checkFrameshifts(const Labeling &,const GffTranscript &,
			const String &substrate,Essex::CompositeNode*);
  bool checkFrameshifts(const Labeling &refLab,
			const Labeling &altLab,
			Essex::CompositeNode *status);
  void translate(GffTranscript &refTrans,GffTranscript &altTrans,
		 String &refProtein,String &altProtein);
  bool hasStartCodon(const String &protein);
  bool hasStopCodon(const String &protein);
  bool hasPTC(const String &protein,int &PTCpos);
  static bool geneIsWellFormed(GffTranscript &,
			       const String &substrate,
			       bool &noStart,bool &noStop,
			       bool &PTC,bool &badSpliceSite,
			       Essex::CompositeNode *status,
			       const SignalSensors &,
			       const int nmdDistParm);
protected:
  const SignalSensors &sensors;
  const GffTranscript &refTrans, &altTrans;
  const Sequence &refSeq, &altSeq;
  const String &refSubstrate, &altSubstrate;
  const Labeling &labeling;
  bool checkDonor(GffExon &refExon,GffExon &altExon,bool &weakened,
		  String &altDonor,String &altWindow,
		  float &refScore,float &altScore,float &cutoff);
  bool checkAcceptor(GffExon &refExon,GffExon &altExon,bool &weakened,
		     String &altDonor,String &altWindow,
		     float &refScore,float &altScore,float &cutoff);
  bool checkDonorStrength(GffExon &refExon,const String &refSubstrate,
		  GffExon &altExon,const String &altSubstrate);
  bool checkAcceptorStrength(GffExon &refExon,const String &refSubstrate,
			     GffExon &altExon,const String &altSubstrate);
  static bool checkDonor(GffExon &exon,const String &substrate,
			 Essex::CompositeNode*,const SignalSensors &);
  static bool checkAcceptor(GffExon &exon,const String &substrate,
			    Essex::CompositeNode*,const SignalSensors &);
  static String getDonor(GffExon &,const String &substrate,int &pos);
  static String getAcceptor(GffExon &,const String &substrate,int &pos);
  float scoreDonor(GffExon &exon,const String &str,const Sequence &seq);
  float scoreAcceptor(GffExon &exon,const String &str,const Sequence &seq);
  String getParsedWindow(SignalSensor &,int windowPos,const String &genome);
  String spliceSiteChangeString(SignalSensor *sensor,
				int refBegin,
				const String &refSubstrate,
				float refScore,
				int altBegin,
				const String &altSubstrate,
				float altScore);
};



