/****************************************************************
 ContentSensor.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ContentSensor_H
#define INCL_ContentSensor_H
#include "BOOM/Sequence.H"
#include "BOOM/String.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/Vector.H"
#include "BOOM/Set.H"
#include "ContentType.H"
#include "BOOM/Strand.H"
class SignalQueue;

class ContentSensor
{
  ContentType contentType; // like INITIAL-EXON, INTRON, 5'-UTR, etc.
  Strand strand;
  BOOM::Set<SignalQueue*> signalQueues;
  static BOOM::Regex binmodRegex;
  
  static ContentSensor *loadBinary(const BOOM::String &);
protected:
  ContentSensor();
  void setContentType(ContentType);
  void setStrand(Strand);
public:
  static const int NO_PHASE;

  // Static methods:
  static ContentSensor *load(const BOOM::String &filename);
  static ContentSensor *load(istream &is);

  // Virtual methods:
  virtual double scoreSingleBase(const Sequence &,const BOOM::String &,
				 int index,Symbol,char)=0;
  virtual void scoreSingleBase(const Sequence &,const BOOM::String &,
			       int index,Symbol,char,double &scorePhase0,
			       double &scorePhase1,double &scorePhase2)=0;
  virtual double scoreSubsequence(const Sequence &,const BOOM::String &,
				  int begin,int length,int phase)=0;
  virtual ContentSensor *reverseComplement()=0;
  virtual bool save(const BOOM::String &filename)=0;
  virtual bool save(ostream &os)=0;
  virtual bool save(BOOM::File &) 
  {throw "binary file format not supported for this model type";}
  virtual void useLogOdds(ContentSensor &nullModel)=0;
  virtual void useLogOdds_anonymous(ContentSensor &nullModel)=0;
  virtual int getOrder()=0;
  virtual bool isPhased()=0;
  virtual ContentSensor *compile() {return this;}
  virtual void reset(const Sequence &,const BOOM::String &,int pos) {}

  // Instance methods:
  virtual ~ContentSensor() {}
  ContentType getContentType() const;
  Strand getStrand() const;
  bool isCoding() const;
  bool isIntron() const;
  bool isIntergenic() const;
  bool isUTR() const;
  void addQueue(SignalQueue *);
  BOOM::Set<SignalQueue*> &getSignalQueues();
  void computeScores(const Sequence &,Array1D<double> &);
  void computeScores(const Sequence &,Array1D<double> &phase0,
		     Array1D<double> &phase1,Array1D<double> &phase2);
};


#endif
