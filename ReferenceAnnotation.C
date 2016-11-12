/****************************************************************
 ReferenceAnnotation.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/GffReader.H"
#include "BOOM/VectorSorter.H"
#include "ReferenceAnnotation.H"
#include "SignalTypeProperties.H"
using namespace std;
using namespace BOOM;

ReferenceAnnotation::ReferenceAnnotation()
  : matrix(NULL), contentRegions(NULL)
{
  // ctor
}



ReferenceAnnotation::ReferenceAnnotation(const String &annoGFFfile,
					 const String &labelingFile,
					 const String &matrixFile,
					 Isochore &isochore,
					 const String &altSeqStr,
					 const Sequence &altSequence)
  : altSeq(altSequence), altSeqStr(altSeqStr)
{
  loadMatrix(matrixFile,isochore.baseDir);
  loadLabeling(labelingFile);
  loadGFF(annoGFFfile,labeling.length());
  //cout<<"CONTENT REGIONS: "<<*contentRegions<<endl;
  initSignals(isochore,altSeqStr,altSequence);
}



ReferenceAnnotation::~ReferenceAnnotation()
{
  delete matrix;
  delete contentRegions;
  // NOTE: don't delete signals here: they are handled by the garbage collector
}



void ReferenceAnnotation::loadMatrix(const String &filename,
				     const String &baseDir)
{
  char *oldPath=new char[PATH_MAX];
  getcwd(oldPath,PATH_MAX);
  chdir(baseDir.c_str());

  matrix=new LabelMatrix(filename);
  matrix->convertToLogs();

  chdir(oldPath);
  delete [] oldPath;
}



const LabelMatrix &ReferenceAnnotation::getMatrix() const
{
  return *matrix;
}



void ReferenceAnnotation::loadLabeling(const String &filename)
{
  labeling.load(filename);
}



const Labeling &ReferenceAnnotation::getLabeling() const
{
  return labeling;
}



const ContentRegions &ReferenceAnnotation::getRegions() const
{
  return *contentRegions;
}



void ReferenceAnnotation::loadGFF(const String &filename,const int seqLen)
{
  Vector<GffTranscript*> &transcripts=*GffReader::loadTranscripts(filename);
  if(transcripts.size()!=1) 
    throw "ReferenceAnnotation::loadGFF() : file must contain exactly one transcript";
  GffTranscript &transcript=*transcripts[0];

  // Add artificial UTR if none is annotated: this is necessary so that we
  // have a TSS and TES signal in the graph
  Vector<GffExon*> utr5, utr3;
  transcript.getUTR(utr5,utr3);
  int begin=transcript.getBegin(), end=transcript.getEnd();
  const int MARGIN=100;
  if(utr5.empty()) 
    transcript.addUTR(new GffExon(ET_SINGLE_UTR5,begin-MARGIN,begin,
				  transcript,false,0,false,0));
  if(utr3.empty()) 
    transcript.addUTR(new GffExon(ET_SINGLE_UTR3,end,end+MARGIN,
				  transcript,false,0,false,0));

  contentRegions=new ContentRegions(transcript,seqLen);
  delete &transcripts;
}



void ReferenceAnnotation::initSignals(Isochore &isochore,
				      const String &altSeqStr,
				      const Sequence &altSequence)
{
  Map<SignalType,SignalSensor*> &sensors=isochore.signalTypeToSensor;
  const Vector<ContentRegion> &regions=contentRegions->asVector();
  const int N=regions.size();
  for(int i=0 ; i<N ; ++i) {
    const ContentRegion &region=regions[i];
    ContentType contentType=region.getType();
    const Interval &interval=region.getInterval();
    const Set<SignalType> &rightSignals=::rightSignals(contentType);
    if(interval.getEnd()<altSeqStr.length())
      if(rightSignals.size()==1) {
	const SignalType rightSignalType=rightSignals.getSingleElement();
	makeSignal(rightSignalType,interval.getEnd(),altSeqStr,
		   altSequence,*sensors[rightSignalType]);
      }
      else switch(contentType) {
	case INTERGENIC:
	case UTR5_SINGLE:
	case UTR5_FINAL:
	  if(i+1<N) {
	    const ContentRegion &nextRegion=regions[i+1];
	    const ContentType nextType=nextRegion.getType();
	    const Set<SignalType> &leftSignals=::leftSignals(nextType);
	    if(leftSignals.size()!=1) INTERNAL_ERROR;
	    const SignalType leftSignalType=leftSignals.getSingleElement();
	    makeSignal(leftSignalType,interval.getEnd(),altSeqStr,
		       altSequence,*sensors[leftSignalType]);
	  }
	  break;
	default: INTERNAL_ERROR;
      }
  }
  sortSignals();
}



SignalPtr ReferenceAnnotation::getStartCodon() const
{
  for(Vector<Signal*>::const_iterator cur=signals.begin(), end=signals.end() ;
      cur!=end ; ++cur)
    if((*cur)->getSignalType()==ATG) return (*cur);
  return NULL;
  //throw "Start codon not found in ReferenceAnnotation::getStartCodon()";
}



void ReferenceAnnotation::sortSignals()
{
  SignalPosComparator signalCmp;
  VectorSorter<Signal*> sorter(signals,signalCmp);
  sorter.sortAscendInPlace();
}



void ReferenceAnnotation::makeSignal(SignalType signalType,int intervalBeginOrEnd,
				     const String &str,const Sequence &seq,
				     SignalSensor &sensor)
{
  int consensusPos, offset=sensor.getConsensusOffset();
  switch(signalType) {
  case ATG: consensusPos=intervalBeginOrEnd;   break;
  case TAG: consensusPos=intervalBeginOrEnd-3; break;
  case GT:  
  case UTR5GT:
  case UTR3GT:
            consensusPos=intervalBeginOrEnd;   break;
  case AG:  
  case UTR5AG:
  case UTR3AG:
            consensusPos=intervalBeginOrEnd-2; break;
  case TSS:
  case TES:
            consensusPos=intervalBeginOrEnd;   break;
  default: INTERNAL_ERROR;
  }
  if(sensor.consensusOccursAt(str,consensusPos))
    makeSignal(signalType,str,seq,consensusPos-offset,sensor);
  else cout<<"WARNING: "<<signalType<<" reference signal does not match: "<<str.substring(consensusPos,2)<<" at "<<consensusPos<<endl; // ### debugging
}



void ReferenceAnnotation::makeSignal(SignalType signalType,const String &str,
				     const Sequence &seq,int contextWindowPos,
				     SignalSensor &sensor)
{
  const double score=sensor.getLogP(seq,str,contextWindowPos);
  if(SignalTypeProperties::global.isEnabled(signalType)) {
    Signal *signal=new Signal(contextWindowPos,score,sensor,sensor.getGC(),
			      signalType);
    signal->setAnnotated();
    signals.push_back(signal);
  }
  if(signalType==GT || signalType==UTR3GT) {
    if(SignalTypeProperties::global.isEnabled(UTR5GT)) {
      Signal *signal=new Signal(contextWindowPos,score,sensor,sensor.getGC(),
				UTR5GT);
      signal->setAnnotated();
      signals.push_back(signal);
    }
    if(SignalTypeProperties::global.isEnabled(UTR3GT)) {
      Signal *signal=new Signal(contextWindowPos,score,sensor,sensor.getGC(),
				UTR3GT);
      signal->setAnnotated();
      signals.push_back(signal);
    }
  }
  else if(signalType==AG || signalType==UTR3AG) {
    if(SignalTypeProperties::global.isEnabled(UTR5AG)) {
      Signal *signal=new Signal(contextWindowPos,score,sensor,sensor.getGC(),
				UTR5AG);
      signal->setAnnotated();
      signals.push_back(signal);
    }
    if(SignalTypeProperties::global.isEnabled(UTR3AG)) {
      Signal *signal=new Signal(contextWindowPos,score,sensor,sensor.getGC(),
				UTR3AG);
      signals.push_back(signal);
    }
  }
}



const Vector<Signal*> &ReferenceAnnotation::getSignals() const
{
  return signals;
}



int ReferenceAnnotation::getStartPosition() const
{
  const Vector<ContentRegion> &regions=contentRegions->asVector();
  for(Vector<ContentRegion>::const_iterator cur=regions.begin(), end=
	regions.end() ; cur!=end ; ++cur) {
    const ContentRegion &region=*cur;
    if(isCoding(region.getType())) return region.getInterval().getBegin();
  }
  throw "Can't find start position in ReferenceAnnotation::getStartPosition()";
}


