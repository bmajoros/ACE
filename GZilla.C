/****************************************************************
 GeneZilla
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <fstream>
#include "GZilla.H"
#include "BOOM/FastaReader.H"
#ifdef REPORT_PROGRESS
#include "BOOM/Progress.H"
#include "BOOM/VectorSorter.H"
#endif

GeneZilla::GeneZilla(const BOOM::String &PROGRAM_NAME,
		     const BOOM::String &VERSION,EdgeFactory &edgeFactory,
		     int &transcriptId)
  : intergenicSums(0),
    seq(NULL),
    seqStr(NULL),
    transcriptId(transcriptId),
    edgeFactory(edgeFactory),
#ifdef EXPLICIT_GRAPHS
    parseGraph(*this),
    oneTerminusOnly(false),
#endif
    PROGRAM_NAME(PROGRAM_NAME),
    VERSION(VERSION),
    isochores(garbageCollector),
    modelCpGislands(false),
    GCregex("gc=(\\d+)"),
    nextIsochoreIndex(0),
    evidenceFilter(NULL),
    useSignalScores(true), useContentScores(true), useDurationScores(true), prohibitPTCs(true) 
  {
    // ctor
    
    ++transcriptId;
    for(int i=0 ; i<3 ; ++i) recentlyEclipsedPhases[i]=false;
  }



GeneZilla::~GeneZilla()
{
  // dtor

  /*
  // Delete the signal sensors
  BOOM::Vector<SignalSensor*>::iterator sscur=signalSensors.begin(), 
    ssend=signalSensors.end();
  for(; sscur!=ssend ; ++sscur)
      delete *sscur;
  */

  // Delete the signal queues
  BOOM::Vector<SignalQueue*>::iterator sqCur=signalQueues.begin(),
    sqEnd=signalQueues.end();
  for(; sqCur!=sqEnd ; ++sqCur) delete *sqCur;

  /*
  // Delete the length distributions
  BOOM::Vector<DiscreteDistribution*>::iterator dCur=distributions.begin(),
    dEnd=distributions.end();
  for(; dCur!=dEnd ; ++dCur)
    delete *dCur;
  */
}



#ifdef EXPLICIT_GRAPHS
ParseGraph &GeneZilla::parse(const BOOM::String &fastaFilename,
			    const BOOM::String &isochoreFilename,
			    float &gcContent)
{
  ParseGraph &graph=parse(fastaFilename,isochoreFilename,seq,seqStr,
			  gcContent);
  delete seq;
  delete seqStr;
  return graph;
}
#endif



#ifdef EXPLICIT_GRAPHS
ParseGraph &GeneZilla::parse(const BOOM::String &fastaFilename,
			    const BOOM::String &isochoreFilename,
			    const Sequence *&substrate,
			    const BOOM::String *&substrateString,
			    float &gcContent)
{
    cerr << "Loading sequence from " << fastaFilename 
	 << "........................" << endl;
    alphabet=DnaAlphabet::global();
    substrate=Sequence::load(fastaFilename,alphabet,substrateId);
    substrateString=substrate->toString(alphabet);
    seqLen=substrate->getLength();

    cerr << "Identifying isochore..." << endl;
    gcContent=getGCcontent(*substrateString);
    cerr << "\t(G+C)/(A+T+C+G)=" << gcContent << endl;

    cerr << "Processing config file..." << endl;
    processIsochoreFile(isochoreFilename,gcContent);

    cerr << "Running the GeneZilla algorithm..." << endl;
    cerr << "\tsequence length=" << substrate->getLength() << endl;
    buildParseGraph(*substrate,*substrateString);

#ifdef REPORT_MEMORY_USAGE
    MemoryProfiler::report("GeneZilla TOTAL MEMORY USAGE: ",cerr);
#endif
    
    return parseGraph;
}
#endif



#ifdef EXPLICIT_GRAPHS
ParseGraph &GeneZilla::parse(const BOOM::String &isochoreFilename,
			    const Sequence &substrate,
			    const BOOM::String &substrateString,
			    float &gcContent)
{
  alphabet=DnaAlphabet::global();
  seqLen=substrate.getLength();
  seq=&substrate;
  seqStr=&substrateString;

  cerr << "Identifying isochore..." << endl;
  gcContent=getGCcontent(substrateString);
  cerr << "\t(G+C)/(A+T+C+G)=" << gcContent << endl;
  
  cerr << "Processing config file..." << endl;
  processIsochoreFile(isochoreFilename,gcContent);
  
  cerr << "Running the GeneZilla algorithm..." << endl;
  cerr << "\tsequence length=" << substrate.getLength() << endl;
  buildParseGraph(substrate,substrateString);
  
#ifdef REPORT_MEMORY_USAGE
  MemoryProfiler::report(
    "                            GeneZilla TOTAL MEMORY USAGE: ",cerr);
#endif
  
  return parseGraph;
}
#endif



#ifdef EXPLICIT_GRAPHS
void GeneZilla::buildParseGraph(const Sequence &seq,const BOOM::String &str)
{
  // Instantiate one signal of each type at the left terminus to act as
  // anchors to which real signals can link back
  instantiateLeftTermini();

  // Make a single left-to-right pass across the sequence
  intergenicSums.resize(seqLen);
  const char *charPtr=str.c_str();
  computeIntergenicSums(seq,str,charPtr);
  for(int pos=0 ; pos<seqLen ; ++pos, ++charPtr)
    {
      Symbol base=seq[pos];
      
      // Check whether any signals occur here
      BOOM::Vector<SignalSensor*>::iterator cur=
	isochore->signalSensors.begin(),
	end=isochore->signalSensors.end();
      for(; cur!=end ; ++cur )
	{
	  SignalSensor &sensor=**cur;
	  if(pos+sensor.getContextWindowLength()>seqLen) continue;

	  SignalPtr signal=sensor.detect(seq,str,pos);
#ifdef FORCE_SPECIFIC_SIGNALS
	  if(!signal && forcedSignalCoords.isMember
	     (pos+sensor.getConsensusOffset()))
	    signal=sensor.detectWithNoCutoff(seq,str,pos);
#endif
	  if(signal) {
	    int begin=signal->getConsensusPosition();
	    int end=signal->posOfBaseFollowingConsensus();
	    bool supported=false;
	    if(evidenceFilter)
	      switch(signal->getSignalType()) 
		{
		case ATG:
		case TAG:
		case NEG_ATG:
		case NEG_TAG:
		  supported=evidenceFilter->codingSignalSupported(begin,end);
		  break;
		case GT:
		case NEG_AG:
		  supported=evidenceFilter->spliceOutSupported(begin);
		  break;
		case AG:
		case NEG_GT:
		  supported=evidenceFilter->spliceInSupported(end);
		  break;
		}
	    else supported=true;

	      if(supported) {
		// Find optimal predecessor for this signal in all 3 phases
		linkBack(str,signal);
		
		// Add this signal to the appropriate queue(s)
		enqueue(signal);
	      }	
	    }
	}
      
      // Check for stop codons & terminate reading frames when they are 
      // found.
      // This check lags behind by two bases so that any stop codon we find 
      // won't overlap with a signal yet to be identified 
      // (consider TAGT=TAG+GT; the TAG should not stop any reading frames 
      // for the GT because when GT is used as a donor only the TA would be 
      // present during translation)
      if(pos>1) handleStopCodons(str,pos-2);

      // Propagate scores of all non-eclipsed signals up to this base
      updateAccumulators(seq,str,pos,base,*charPtr);
    }

  // Instantiate an anchor signal of each type at the right terminus
  // and link them back, to complete the dynamic programming evaluation:
  double parseScore;
  BOOM::Stack<SignalPtr> *path=
    instantiateRightTermini(str,seqLen,parseScore);
  delete path;

  // Run the garbage collector & build the parse graph
  BOOM::Vector<SignalPtr>::iterator rCur=rightTermini.begin(), 
    rEnd=rightTermini.end();
  for(; rCur!=rEnd ; ++rCur)
    {
      SignalPtr s=*rCur;
      garbageCollector.markLeft(s);
    }
  BOOM::Vector<SignalPtr>::iterator lCur=leftTermini.begin(), 
    lEnd=leftTermini.end();
  for(; lCur!=lEnd ; ++lCur)
    {
      SignalPtr s=*lCur;
      garbageCollector.markRight(s);
    }
  garbageCollector.sweep();
  BOOM::Set<SignalPtr>::iterator sCur=garbageCollector.signalsBegin(),
    sEnd=garbageCollector.signalsEnd();
  while(sCur!=sEnd)
    {
      SignalPtr s=*sCur;
      if(!useSignalScores) s->dropSignalScores();
      if(!useContentScores) s->dropContentScores();
      parseGraph.addVertex(s);
      ++sCur;
      garbageCollector.drop(s);
    }
}
#endif



int GeneZilla::main(int argc,char *argv[])
  {
    // Process command line
    BOOM::CommandLine cmd(argc,argv,"r:");
    if(cmd.numArgs()!=2)
      throw string(
"\ngenezilla <*.iso> <*.fasta> [-r <*.gff>]\n\
          where -r <*.gff> specifies a GFF file to load and score\n\n");
    BOOM::String isochoreFilename=cmd.arg(0);
    BOOM::String fastaFilename=cmd.arg(1);

    cerr << "Loading sequence from " << fastaFilename 
	 << "........................" << endl;
    alphabet=DnaAlphabet::global();

    BOOM::FastaReader fastaReader(fastaFilename);
    BOOM::String defline, substrateString;
    fastaReader.nextSequence(defline,substrateString);
    Sequence substrate(substrateString,DnaAlphabet::global());

    seqLen=substrate.getLength();

    /*
    if(cmd.option('r')) 
      loadGFF(cmd.optParm('r'),substrate,substrateString);
    */

    cerr << "Running the GeneZilla algorithm..." << endl;
    cerr << "\tsequence length=" << substrate.getLength() << endl;
    ofstream osGraph;
    BOOM::Stack<SignalPtr> *path=mainAlgorithm(substrate,substrateString,
					       osGraph,false,"");
    delete path;

#ifdef REPORT_MEMORY_USAGE
    MemoryProfiler::report("TOTAL MEMORY USAGE: ",cerr);
#endif

    return 0;
  }



BOOM::Stack<SignalPtr> * GeneZilla::processChunk(const Sequence &substrate,
					const BOOM::String &substrateString,
					const BOOM::String &isoFilename,
					const BOOM::String &substrateId,
						 ostream &osGraph,
						 bool dumpGraph,
						 String psaFilename)
{
  seq=&substrate;
  seqStr=&substrateString;
  seqLen=substrate.getLength();
  this->substrateId=substrateId;
  
  nextIsochoreInterval.begin=-1;
  if(!isochoreIntervals.isDefined(substrateId))
    gcContent=getGCcontent(substrateString);
  else
    {
      BOOM::Vector<IsochoreInterval> &intervals=
	isochoreIntervals[substrateId];
      if(intervals.size()>0)
	{
	  IsochoreInterval &interval=intervals[0];
	  gcContent=interval.GC;
	  nextIsochoreIndex=1;
	  if(intervals.size()>1)
	    nextIsochoreInterval=intervals[1];
	  else
	    nextIsochoreInterval.begin=-1;
	}
      else
	gcContent=getGCcontent(substrateString);
    }

  cerr<<"GC content = "<<gcContent<<endl;
  if(isochores.getNumIsochores()==0)
    {
      // This is the first chunk we are seeing...so do some initialization
      // first:
      cerr << "Processing config file..." << endl;
      processIsochoreFile(isoFilename,gcContent);
    }
  else 
    {
      // This is not the first chunk...just switch the isochore:
      switchIsochore(gcContent,0);
      resetAccumulatorPositions();
      BOOM::Vector<SignalQueue*>::iterator qCur=signalQueues.begin(),
	qEnd=signalQueues.end();
      for(; qCur!=qEnd ; ++qCur)
	{
	  SignalQueue *queue=*qCur;
	  queue->resetQueue(isochore);
	}
    }

  cerr << "\tsequence length=" << seqLen << endl;
  return mainAlgorithm(substrate,substrateString,osGraph,dumpGraph,
		       psaFilename);
}



float GeneZilla::getGCcontent(const BOOM::String &seq)
{
  int n=seq.length(), ATCG=0, GC=0;
  const char *p=seq.c_str();
  for(int i=0 ; i<n ; ++i)
    switch(*p++)
      {
      case 'G':
      case 'C':
	++GC;
	// no break...fall through:
      case 'A':
      case 'T':
	++ATCG;
      }
  return GC/float(ATCG);
}



void GeneZilla::instantiateLeftTermini()
{
  SignalSensor &sensor=*isochore->signalTypeToSensor[LEFT_TERMINUS];
  SignalPtr terminus=sensor.getLeftTerminus(0);
  Set<ContentType> &queues=sensor.belongsInWhichQueues();
  for(Set<ContentType>::iterator cur=queues.begin(), end=queues.end() ;
      cur!=end ; ++cur) {
    ContentType contentType=*cur;
    SignalQueue *queue=contentToQueue[contentType];
    queue->addSignal(terminus);
    for(int i=0 ; i<3 ;++i) terminus->getInductiveScore(i)=0;
  }
#ifdef EXPLICIT_GRAPHS
  leftTermini.push_back(terminus);
#endif
}



#ifdef EXPLICIT_GRAPHS
void GeneZilla::useOneTerminusOnly()
{
  oneTerminusOnly=true;
}
#endif



BOOM::Stack<SignalPtr> *GeneZilla::instantiateRightTermini(
						    const BOOM::String &str,
						    int seqLen,
						    double &bestScore)
{
  int bestPhase=INT_MIN;
  bestScore=NEGATIVE_INFINITY;
  SignalSensor &sensor=*isochore->signalTypeToSensor[RIGHT_TERMINUS];
  BOOM::Set<ContentType> &queues=sensor.linksBackToWhichQueues();
  ContentType contentType=*queues.begin();
  SignalPtr rightTerminus=sensor.getRightTerminus(seqLen,0);
  linkBack(str,rightTerminus);
  Propagator &prop=rightTerminus->getPropagator(contentType);
  for(int i=0 ; i<3 ; ++i) {
    double thisScore=prop[i];
    if(thisScore>bestScore) {
      bestScore=thisScore;
      bestPhase=i;
    }
  }
#ifdef EXPLICIT_GRAPHS
  rightTermini.push_back(rightTerminus);
#endif
  return traceBack(rightTerminus,bestPhase);
}



BOOM::Stack<SignalPtr> * GeneZilla::mainAlgorithm(const Sequence &seq,
						  const BOOM::String &str,
						  ostream &osGraph,
						  bool dumpGraph,
						  String psaFilename)
{
  throw "not used";

  /*
  // Compute cumulative intergenic score at each base
  const char *charPtr=str.c_str();
  intergenicSums.resize(seqLen);
  computeIntergenicSums(seq,str,charPtr);
  if(psaFilename.length()>0) {
    ofstream os(psaFilename.c_str());
    for(int i=0 ; i<seqLen ; ++i)
      os<<intergenicSums[i]<<endl;
  }

#ifndef EXPLICIT_GRAPHS
  // Instantiate one signal of each type at the left terminus to act as
  // anchors to which real signals can link back
  instantiateLeftTermini();

#ifdef REPORT_PROGRESS
  BOOM::Progress progress;
  progress.start(seqLen);
#endif

  // Make a single left-to-right pass across the sequence
  for(int pos=0 ; pos<seqLen ; ++pos, ++charPtr)
    {
      if(pos==nextIsochoreInterval.begin) crossIsochoreBoundary(pos);
      Symbol base=seq[pos];
      
      // Check whether any signals occur here
      BOOM::Vector<SignalSensor*>::iterator cur=
	isochore->signalSensors.begin(), end=isochore->signalSensors.end();
      for(; cur!=end ; ++cur )
	{
	  SignalSensor &sensor=**cur;
	  if(pos+sensor.getContextWindowLength()>seqLen) continue;
#ifdef FORCE_SPECIFIC_SIGNALS
	  SignalPtr signal=
	    (forcedSignalCoords.isMember(pos+sensor.getConsensusOffset()) ? 
	     sensor.detectWithNoCutoff(seq,str,pos) :
	     sensor.detect(seq,str,pos));
#else
	  SignalPtr signal=sensor.detect(seq,str,pos);
#endif
	  if(signal)
	    {
	      // Find optimal predecessor for this signal in all 3 phases
	      linkBack(str,signal);

	      // Add this signal to the appropriate queue(s)
	      enqueue(signal);
	    }
	}

      // Check for stop codons & terminate reading frames when they 
      // are found.  This check lags behind by two bases so that any 
      // stop codon we find won't overlap with a signal yet to be 
      // identified (consider TAGT=TAG+GT; the TAG should not stop 
      // any reading frames for the GT because when GT is used as a 
      // donor only the TA would be present during translation)
      if(pos>1) handleStopCodons(str,pos-2);

      // Propagate scores of all non-eclipsed signals up to this base
      updateAccumulators(seq,str,pos,base,*charPtr);

#ifdef REPORT_PROGRESS
      if((pos+1)%250000==0) cerr<<progress.getProgress(pos)<<endl;
#endif
    }

  // Instantiate an anchor signal of each type at the right terminus
  // and link them back, to complete the dynamic programming evaluation:
  double parseScore;
  BOOM::Stack<SignalPtr> *path=
    instantiateRightTermini(str,seqLen,parseScore);

  // Output gene prediction in GFF format
  generateGff(path,seqLen,parseScore);
#endif

#ifdef EXPLICIT_GRAPHS
  //###DEBUGGING: do the prediction using the graph instead of the "trellis"
  //const char *charPtr=str.c_str();
  //intergenicSums.resize(seqLen);
  //computeIntergenicSums(seq,str,charPtr);
  buildParseGraph(seq,str);
  double parseScore;
  BOOM::Stack<SignalPtr> *path=parseGraph.findOptimalPath(parseScore);
  generateGff(path,seqLen,parseScore);
  if(dumpGraph) {
    parseGraph.setVertexIndices();
    osGraph<<parseGraph<<endl;
  }
#endif

#ifdef REPORT_MEMORY_USAGE
    MemoryProfiler::report("GeneZilla TOTAL MEMORY USAGE: ",cerr);
#endif

  return path;
  */
}



void GeneZilla::updateAccumulators(const Sequence &seq,
				  const BOOM::String &str,
				  int pos,Symbol base,char c)
{
  double score, scorePhase0, scorePhase1, scorePhase2;
  BOOM::Vector<ContentSensor*>::iterator cur=
    isochore->contentSensors.begin(),
    end=isochore->contentSensors.end();
  for(; cur!=end ; ++cur)
    {
      ContentSensor &contentSensor=**cur;
      bool isCoding=contentSensor.isCoding();
      if(isCoding)
	contentSensor.scoreSingleBase(seq,str,pos,base,c,scorePhase0,
				      scorePhase1,scorePhase2);
      else
	score=contentSensor.scoreSingleBase(seq,str,pos,base,c);
      BOOM::Set<SignalQueue*> &queues=contentSensor.getSignalQueues();
      BOOM::Set<SignalQueue*>::iterator cur=queues.begin(), 
	end=queues.end();
      for(; cur!= end ; ++cur)
	{
	  SignalQueue &queue=**cur;
	  if(isCoding)
	    queue.addToAccumulator(scorePhase0,scorePhase1,scorePhase2,pos);
	  else
	    queue.addToAccumulator(score);
	}
    }
}



void GeneZilla::linkBack(const BOOM::String &str,SignalPtr newSignal)
{
  SignalPtr bestPred[3];
  bestPred[0]=bestPred[1]=bestPred[2]=NULL;
  double bestScore[3];
  bestScore[0]=bestScore[1]=bestScore[2]=NEGATIVE_INFINITY;
  int newConsPos=newSignal->getConsensusPosition();
  Strand strand=newSignal->getStrand();
  SignalType signalType=newSignal->getSignalType();
  if(::endsCoding(signalType)) observeRecentStopCodons(str,newSignal);

  // Consider all queues that newSignal could link back to
  BOOM::Set<ContentType> &queues=newSignal->linksBackToWhichQueues();
  BOOM::Set<ContentType>::iterator cur=queues.begin(), end=queues.end();
  for(; cur!=end ; ++cur)
    {
      // Make sure the queue's signals are all propagated up to this point
      ContentType contentType=*cur;
      SignalQueue *queue=contentToQueue[contentType];
      queue->flushAccumulator();

      // Consider all of the queue's signals as potential predecessors
      selectPredecessors(newConsPos,*queue,contentType,strand,
			 bestScore,bestPred,signalType,str,newSignal);
    }

  // Install the selected predecessors as dyn. prog. links and update
  // inductive scores for all propagators of this signal (identically)
  BOOM::Set<ContentType> &nextQueues=newSignal->belongsInWhichQueues();
  cur=nextQueues.begin();
  end=nextQueues.end();
  for(; cur!=end ; ++cur)
    {
      ContentType contentType=*cur;
      Propagator &prop=newSignal->getPropagator(contentType);
      prop[0]+=bestScore[0];
      prop[1]+=bestScore[1];
      prop[2]+=bestScore[2];
      
      // At this point the propagators are all the same, so if we overwrite
      // one set of inductive scores with those for another content type, it
      // won't matter, because they are the same:
      newSignal->getInductiveScore(0)=prop[0];
      newSignal->getInductiveScore(1)=prop[1];
      newSignal->getInductiveScore(2)=prop[2];
    }
  newSignal->setPredecessor(0,bestPred[0]);
  newSignal->setPredecessor(1,bestPred[1]);
  newSignal->setPredecessor(2,bestPred[2]);
}



inline void GeneZilla::selectPredecessors(int newConsPos,
					 SignalQueue &queue,
					 ContentType contentType,
					 Strand strand,
					 double bestScore[3],
					 SignalPtr bestPred[3],
					 SignalType toType,
					 const BOOM::String &substrate,
					 SignalPtr signal)
{
  switch(contentType)
    {
    case UTR5_INITIAL:
    case UTR5_INTERNAL:
    case UTR5_FINAL:
    case UTR5_SINGLE:
    case NEG_UTR5_INITIAL:
    case NEG_UTR5_INTERNAL:
    case NEG_UTR5_FINAL:
    case NEG_UTR5_SINGLE:
    case UTR3_INITIAL:
    case UTR3_INTERNAL:
    case UTR3_FINAL:
    case UTR3_SINGLE:
    case NEG_UTR3_INITIAL:
    case NEG_UTR3_INTERNAL:
    case NEG_UTR3_FINAL:
    case NEG_UTR3_SINGLE:
    case INTERGENIC:
      selectIntergenicPred(newConsPos,queue,strand,bestScore,bestPred,
			   contentType,toType,signal);
      break;
    case INTRON:
    case UTR5_INTRON:
    case UTR3_INTRON:
    case NEG_INTRON:
    case NEG_UTR5_INTRON:
    case NEG_UTR3_INTRON:
      selectIntronPred(newConsPos,queue,strand,bestScore,bestPred,
		       contentType,toType,substrate,signal);
      break;
    default:
      selectCodingPred(newConsPos,queue,strand,bestScore,bestPred,
		       contentType,toType,signal);
    }
}



void GeneZilla::selectCodingPred(int newConsPos,
				SignalQueue &queue,
				Strand strand,double bestScore[3],
				SignalPtr bestPred[3],
				ContentType contentType,
				SignalType toType,
				SignalPtr signal)
{
  // Version #1 of 3: CODING.  Phase is translated across the exon.
  BOOM::Iterator<SignalPtr> &cur=queue.begin(), &end=queue.end();
  for(; cur!=end ; ++cur) {
    SignalPtr pred=*cur;
    Propagator &predProp=pred->getPropagator(contentType);
    SignalType predType=pred->getSignalType();
    int oldPos=pred->posOfBaseFollowingConsensus();
    int length=newConsPos-oldPos;
    if(length<0) continue;
    double lengthScore=useDurationScores ? 
      isochore->contentToDistribution[contentType]->
      getLogP(length+pred->getConsensusLength()) : 0;
    int frameDelta=length % 3;
#ifdef EXPLICIT_GRAPHS
    double linkScores[3];
#endif
    for(int oldPhase=0 ; oldPhase<3 ; ++oldPhase) {
      if(isinf(predProp[oldPhase])) {
#ifdef EXPLICIT_GRAPHS
	linkScores[oldPhase]=NEGATIVE_INFINITY;
#endif
	continue;
      }
      int newPhase=(strand==FORWARD_STRAND ?
		    (oldPhase+frameDelta) % 3 :
		    posmod(oldPhase-frameDelta));
      if(recentlyEclipsedPhases[newPhase]) {
#ifdef EXPLICIT_GRAPHS
	if(prohibitPTCs) linkScores[oldPhase]=NEGATIVE_INFINITY;
#endif
	continue;
      }

      // Implement phase-specific intron states:
      double transScore=
	scoreIntronPhases(predType,toType,oldPhase,newPhase);

      double predScore=
	predProp[oldPhase] + lengthScore + transScore;

      double &bestPredScore=bestScore[newPhase];

#ifdef EXPLICIT_GRAPHS
      linkScores[oldPhase]=
	predScore-pred->posteriorInductiveScore(oldPhase);
#endif
      if(finite(predScore) && predScore>bestPredScore) {
	bestPredScore=predScore;
	bestPred[newPhase]=pred;
      }
    }
#ifdef EXPLICIT_GRAPHS

    if(finite(linkScores[0]) || finite(linkScores[1]) || 
       finite(linkScores[2])) {
      // ### NEED TO CHECK PHASE IF TAG/ATG INVOLVED
      SignalType t=signal->getSignalType();
      if(t==ATG || t==NEG_TAG) {
	int predPhase;
	switch(strand) {
	case FORWARD_STRAND:
	  predPhase=posmod(0-length);
	case REVERSE_STRAND:
	  predPhase=(0+length)%3;
	}
	if(!isFinite(linkScores[predPhase])) continue;
	linkScores[(predPhase+1)%3]=linkScores[(predPhase+2)%3]=
	  NEGATIVE_INFINITY;
      }
      
      Edge *edge=
	edgeFactory.newPhasedEdge(linkScores[0],linkScores[1],
				  linkScores[2],pred,signal);
      if(edge) {
	pred->addEdgeOut(edge);
	signal->addEdgeIn(edge);
      }	
    }
#endif
  }
}



inline void GeneZilla::selectIntronPred(int newConsPos,SignalQueue &q,
				       Strand strand,double bestScore[3],
				       SignalPtr bestPred[3],
				       ContentType contentType,
				       SignalType toType,
				       const BOOM::String &substrate,
				       SignalPtr signal)
{
  // Version #2 of 3: INTRONS.  Phase is not translated across the intron.
  // Also, we look for stop codons interrupted by the intron and close their
  // open reading frames.

  SignalSensor *stopCodonSensor=isochore->stopCodonSensor;
  SignalSensor *negStopCodonSensor=isochore->negStopCodonSensor;
  int seqLen=substrate.length();
  bool phaseIsOpen[3];
  for(int i=0 ; i<3 ; ++i) phaseIsOpen[i]=true;
  BOOM::String nextCodingBase,nextTwoCodingBases;
  if(newConsPos<seqLen-2) {
    // ### this would be faster if I didn't use STL strings....
    nextCodingBase=substrate.substring(newConsPos+2,1);
    nextTwoCodingBases=substrate.substring(newConsPos+2,2);
  }
  IntronQueue &queue=(IntronQueue&)q;
  BOOM::Map<SignalPtr,int>::iterator cur=queue.beginUniq(), 
    end=queue.endUniq();
  for(; cur!=end ; ++cur) {
    SignalPtr pred=(*cur).first;
    int oldConsPos=pred->getConsensusPosition();

    // Handle stop codons interrupted by the intron
    if(prohibitPTCs)
      if(newConsPos<seqLen-2 && oldConsPos>1) {
	BOOM::String prevCodingBase=substrate.substring(oldConsPos-1,1);
	BOOM::String prevTwoCodingBases=
	  substrate.substring(oldConsPos-2,2);
	BOOM::String codon12=prevCodingBase+nextTwoCodingBases;
	BOOM::String codon21=prevTwoCodingBases+nextCodingBase;
	phaseIsOpen[0]=
	  strand==FORWARD_STRAND ?
	  true :
	  !negStopCodonSensor->consensusOccursAt(codon21,0);
	phaseIsOpen[1]=
	  strand==FORWARD_STRAND ?
	  !stopCodonSensor->consensusOccursAt(codon12,0) :
	  !negStopCodonSensor->consensusOccursAt(codon12,0);
	phaseIsOpen[2]=
	  strand==FORWARD_STRAND ?
	  !stopCodonSensor->consensusOccursAt(codon21,0) :
	  true;
      }

    double transScore=
	isochore->transitionProbs->getLogP(pred->getSignalType(),toType);
    int oldPos=pred->posOfBaseFollowingConsensus();
    int length=newConsPos-oldPos;
    if(length<0) continue;
    double lengthScore=useDurationScores ? 
      isochore->contentToDistribution[contentType]->
      getLogP(length+pred->getConsensusLength()) : 0;
    Propagator &predProp=pred->getPropagator(contentType);
#ifdef EXPLICIT_GRAPHS
    double linkScores[3];
#endif
    for(int phase=0 ; phase<3 ; ++phase) {
      if(!phaseIsOpen[phase]) {
#ifdef EXPLICIT_GRAPHS
	linkScores[phase]=NEGATIVE_INFINITY;
#endif
	continue;
      }
      double predScore=predProp[phase] + lengthScore + transScore;
#ifdef EXPLICIT_GRAPHS
      linkScores[phase]=predScore-pred->posteriorInductiveScore(phase);
#endif
      double &bestPredScore=bestScore[phase];
      if(finite(predScore) && predScore>bestPredScore){
	bestPredScore=predScore;
	bestPred[phase]=pred;
      }
    }
#ifdef EXPLICIT_GRAPHS
      if(finite(linkScores[0]) || finite(linkScores[1]) || 
	 finite(linkScores[2]))	{
	Edge *edge=
	  edgeFactory.newPhasedEdge(linkScores[0],linkScores[1],
				    linkScores[2],pred,signal);
	if(edge) {
	  pred->addEdgeOut(edge);
	  signal->addEdgeIn(edge);
	}
      }
#endif
  }
}



inline void GeneZilla::selectIntergenicPred(int newConsPos,
					   SignalQueue &queue,
					   Strand strand,
					   double bestScore[3],
					   SignalPtr bestPred[3],
					   ContentType contentType,
					   SignalType toType,
					   SignalPtr signal)
{
  // Version #3 of 3: INTERGENIC.  Phase 0 or 2 is used depending on strand.

  int newPhase=(strand==FORWARD_STRAND ? 0 : 2);
  BOOM::Iterator<SignalPtr> &cur=queue.begin(), &end=queue.end();
  for(; cur!=end ; ++cur) {
    SignalPtr pred=*cur;
    double transScore=
      isochore->transitionProbs->getLogP(pred->getSignalType(),toType);
    int oldPos=pred->posOfBaseFollowingConsensus();
    int length=newConsPos-oldPos;
    if(length<0) continue;
    double lengthScore=useDurationScores ? 
      isochore->contentToDistribution[contentType]->
      getLogP(length+pred->getConsensusLength()) : 0;
    Propagator &predProp=pred->getPropagator(contentType);
    int oldPhase=(pred->getStrand()==FORWARD_STRAND ? 0 : 2);
    double predScore=predProp[oldPhase] + lengthScore + transScore;
    if(isUTR3(contentType)) 
      predScore+=isochore->threePrimeOptimism;
    else if(isUTR5(contentType))
      predScore+=isochore->fivePrimeOptimism; 
    double &bestPredScore=bestScore[newPhase];
    if(finite(predScore)) {
#ifdef EXPLICIT_GRAPHS
      Edge *edge=
	edgeFactory.newNonPhasedEdge(
	      predScore-pred->posteriorInductiveScore(oldPhase),pred,
	      signal);
      pred->addEdgeOut(edge);
      signal->addEdgeIn(edge);
#endif
      if(predScore>bestPredScore) {
	bestPredScore=predScore;
	bestPred[newPhase]=pred;
      }
    }
  }
}



void GeneZilla::enqueue(SignalPtr signal)
{
  BOOM::Set<ContentType> &queues=signal->belongsInWhichQueues();
  BOOM::Set<ContentType>::iterator cur=queues.begin(), end=queues.end();
  for(; cur!=end ; ++cur) {
    ContentType contentType=*cur;
    contentToQueue[contentType]->addSignal(signal);
  }
}



void GeneZilla::handleStopCodons(const BOOM::String &str,int pos)
{
  if(!prohibitPTCs) return;
  if(isochore->stopCodonSensor->consensusOccursAt(str,pos))
    terminateForwardORFs(pos);
#ifndef FORWARD_STRAND_ONLY
  else if(isochore->negStopCodonSensor->consensusOccursAt(str,pos))
    terminateReverseORFs(pos);
#endif
}



void GeneZilla::terminateForwardORFs(int TAGposition)
{
  if(!prohibitPTCs) return;

  INTERNAL_ERROR;

  BOOM::Vector<SignalQueue*>::iterator qCur=forwardCodingQueues.begin(),
    qEnd=forwardCodingQueues.end();
  for(; qCur!=qEnd ; ++qCur)
    {
      SignalQueue &queue=**qCur;

      // *********** ITERATE THROUGH THE MAIN QUEUE *************

      BOOM::Iterator<SignalPtr> &sCur=queue.begin(), &sEnd=queue.end();
      while(sCur!=sEnd)
	{
	  // Figure out which phase the TAG occurs in, and terminate it
	  SignalPtr signal=*sCur;
	  int delta=TAGposition-signal->posOfBaseFollowingConsensus();
	  if(delta<0) {++sCur;continue;}
	  int stoppedPhase=posmod(-delta);
	  Propagator &prop=signal->getPropagator(queue.getContentType());

	  prop[stoppedPhase]=NEGATIVE_INFINITY;

	  // If all phases of this signal are stopped, remove it from queue
	  if(isinf(prop[(stoppedPhase+1)%3]) &&
	     isinf(prop[(stoppedPhase+2)%3]))
	    {
	      BOOM::Iterator<SignalPtr> &victim=sCur.clone();
	      ++sCur; // otherwise queue.erase() will invalidate it
	      queue.drop(victim);
	      delete &victim;
	      continue;
	    }
	  ++sCur;
	}

      // ********** NOW DO THE SAME THING FOR THE HOLDING QUEUE *********

      BOOM::List<SignalPtr>::iterator hsCur=queue.getHoldingQueue();
      BOOM::List<SignalPtr>::iterator hsEnd=queue.holdingQueueEnd();
      while(hsCur!=hsEnd)
	{
	  // Figure out which phase the TAG occurs in, and terminate it
	  SignalPtr signal=*hsCur;
	  int delta=TAGposition-signal->posOfBaseFollowingConsensus();
	  if(delta<0) {++hsCur;continue;}
	  int stoppedPhase=posmod(-delta);
	  Propagator &prop=signal->getPropagator(queue.getContentType());

	  prop[stoppedPhase]=NEGATIVE_INFINITY;

	  // If all phases of this signal are stopped, remove it from queue
	  if(isinf(prop[(stoppedPhase+1)%3]) &&
	     isinf(prop[(stoppedPhase+2)%3]))
	    {
	      BOOM::List<SignalPtr>::iterator victim=hsCur;
	      ++hsCur; // otherwise queue.erase() will invalidate it
	      queue.dropFromHoldingQueue(victim);
	      continue;
	    }
	  ++hsCur;
	}
    }
}



void GeneZilla::terminateReverseORFs(int TAGposition)
{
  if(!prohibitPTCs) return;
  BOOM::Vector<SignalQueue*>::iterator qCur=reverseCodingQueues.begin(),
    qEnd=reverseCodingQueues.end();
  for(; qCur!=qEnd ; ++qCur)
    {
      SignalQueue &queue=**qCur;

      // *********** ITERATE THROUGH THE MAIN QUEUE *************

      BOOM::Iterator<SignalPtr> &sCur=queue.begin(), &sEnd=queue.end();
      while(sCur!=sEnd)
	{
	  // Figure out which phase the TAG occurs in, and terminate it
	  SignalPtr signal=*sCur;
	  int delta=TAGposition-signal->posOfBaseFollowingConsensus();
	  if(delta<0) {++sCur;continue;}
	  int stoppedPhase=posmod(delta-1);
	  Propagator &prop=signal->getPropagator(queue.getContentType());
	  prop[stoppedPhase]=NEGATIVE_INFINITY;

	  // If all phases of this signal are stopped, remove it from queue
	  if(isinf(prop[(stoppedPhase+1)%3]) &&
	     isinf(prop[(stoppedPhase+2)%3]))
	    {
	      BOOM::Iterator<SignalPtr> &victim=sCur.clone();
	      ++sCur; // otherwise queue.erase() will invalidate it
	      queue.drop(victim);
	      delete &victim;
	      continue;
	    }
	  ++sCur;
	}

      // ********** NOW DO THE SAME THING FOR THE HOLDING QUEUE *********

      BOOM::List<SignalPtr>::iterator hsCur=queue.getHoldingQueue();
      BOOM::List<SignalPtr>::iterator hsEnd=queue.holdingQueueEnd();
      while(hsCur!=hsEnd)
	{
	  // Figure out which phase the TAG occurs in, and terminate it
	  SignalPtr signal=*hsCur;
	  int delta=TAGposition-signal->posOfBaseFollowingConsensus();
	  if(delta<0) {++hsCur;continue;}
	  int stoppedPhase=posmod(delta-1);
	  Propagator &prop=signal->getPropagator(queue.getContentType());
	  prop[stoppedPhase]=NEGATIVE_INFINITY;

	  // If all phases of this signal are stopped, remove it from queue
	  if(isinf(prop[(stoppedPhase+1)%3]) &&
	     isinf(prop[(stoppedPhase+2)%3]))
	    {
	      BOOM::List<SignalPtr>::iterator victim=hsCur;
	      ++hsCur; // otherwise queue.erase() will invalidate it
	      queue.dropFromHoldingQueue(victim);
	      continue;
	    }
	  ++hsCur;
	}
    }  
}



BOOM::Stack<SignalPtr> *GeneZilla::traceBack(SignalPtr rightTerminus,
					     int phase)
{
  cerr << "Tracing back to find optimal path..." << endl;

  // Trace back and build stack of visited signals
  BOOM::Stack<SignalPtr> *stk=new BOOM::Stack<SignalPtr>;
  if(!rightTerminus) {cout<<"RT=null"<<endl; return stk;} // ###
  stk->push(rightTerminus);
  SignalPtr signal=rightTerminus;
  while(true)
    {
      SignalPtr pred=signal->getPredecessor(phase);
      signal->setPredecessor((phase+1)%3,NULL);
      signal->setPredecessor((phase+2)%3,NULL);
      if(!pred) break;
      int distance=signal->getConsensusPosition() -
	pred->posOfBaseFollowingConsensus();
      stk->push(pred);
      switch(pred->getSignalType())
	{
	case NEG_TAG:
	case NEG_TES:
	case NEG_ATG:     
	case NEG_TSS:    
	  phase=2;
	  break;
	case TSS:
	case ATG:
	case TAG:         
	case TES:       
	  phase=0;
	  break;
	case GT:
	case NEG_AG:
	  break;
	case AG:
	  phase=posmod(phase-distance);
	  break;
	case NEG_GT:
	  phase=posmod(phase+distance);
	  break;
	}
      signal=pred;
    }
  return stk;
}



int GeneZilla::getPhase(SignalPtr signal)
{
  // ### This needs to be fixed to work for partial genes

  int phase=(signal->getStrand()==FORWARD_STRAND ? 0 : 2);
  for(int i=0 ; i<3 ; ++i)
    if(signal->getPredecessor(i))
      phase=i;
  return phase;
}



double GeneZilla::scoreIntronPhases(SignalType predType,SignalType toType,
				       int oldPhase,int newPhase)
{
  Transitions *transitionProbs=isochore->transitionProbs;
  double transScore;
  switch(toType)
    {
    case NEG_AG:
      switch(predType)
	{
	case NEG_TAG://                      * * * * * Final exon
	  switch(newPhase)
	    {
	    case 0:
	      transScore=transitionProbs->getLogP(NEG_TAG,NEG_AG0);
	      break;
	    case 1:
	      transScore=transitionProbs->getLogP(NEG_TAG,NEG_AG1);
	      break;
	    case 2:
	      transScore=transitionProbs->getLogP(NEG_TAG,NEG_AG2);
	      break;
	    }
	  break;
	case NEG_GT://                    * * * * * Internal exon
	  switch(newPhase)
	    {
	    case 0:
	      switch(oldPhase)
		{
		case 0: 
		  transScore=transitionProbs->getLogP(NEG_GT0,NEG_AG0);
		  break;
		case 1: 
		  transScore=transitionProbs->getLogP(NEG_GT1,NEG_AG0);
		  break;
		case 2: 
		  transScore=transitionProbs->getLogP(NEG_GT2,NEG_AG0);
		  break;
		}
	      break;
	    case 1:
	      switch(oldPhase)
		{
		case 0: 
		  transScore=transitionProbs->getLogP(NEG_GT0,NEG_AG1);
		  break;
		case 1: 
		  transScore=transitionProbs->getLogP(NEG_GT1,NEG_AG1);
		  break;
		case 2: 
		  transScore=transitionProbs->getLogP(NEG_GT2,NEG_AG1);
		  break;
		}
	      break;
	    case 2:
	      switch(oldPhase)
		{
		case 0: 
		  transScore=transitionProbs->getLogP(NEG_GT0,NEG_AG2);
		  break;
		case 1: 
		  transScore=transitionProbs->getLogP(NEG_GT1,NEG_AG2);
		  break;
		case 2: 
		  transScore=transitionProbs->getLogP(NEG_GT2,NEG_AG2);
		  break;
		}
	      break;
	    }
	  break;
	}
      break;
    case NEG_ATG:
      switch(predType)
	{
	case NEG_TAG://                     * * * * * Single exon
	  transScore=transitionProbs->getLogP(NEG_TAG,NEG_ATG);break;
	  break;
	case NEG_GT://                  * * * * * * Initial exon
	  switch(oldPhase)
	    {
	    case 0:
	      transScore=transitionProbs->getLogP(NEG_GT0,NEG_ATG);break;
	    case 1:
	      transScore=transitionProbs->getLogP(NEG_GT1,NEG_ATG);break;
	    case 2:
	      transScore=transitionProbs->getLogP(NEG_GT2,NEG_ATG);break;
	    }
	  break;
	}
      break;
    case GT:
      switch(predType)
	{
	case ATG://                  * * * * * * Initial exon
	  switch(newPhase)
	    {
	    case 0:
	      transScore=transitionProbs->getLogP(ATG,GT0);break;
	    case 1:
	      transScore=transitionProbs->getLogP(ATG,GT1);break;
	    case 2:
	      transScore=transitionProbs->getLogP(ATG,GT2);break;
	    }
	  break;
	case AG://                    * * * * * Internal exon
	  switch(oldPhase)
	    {
	    case 0:
	      switch(newPhase)
		{
		case 0:
		  transScore=transitionProbs->getLogP(AG0,GT0);break;
		case 1:
		  transScore=transitionProbs->getLogP(AG0,GT1);break;
		case 2:
		  transScore=transitionProbs->getLogP(AG0,GT2);break;
		}
	      break;
	    case 1:
	      switch(newPhase)
		{
		case 0:
		  transScore=transitionProbs->getLogP(AG1,GT0);break;
		case 1:
		  transScore=transitionProbs->getLogP(AG1,GT1);break;
		case 2:
		  transScore=transitionProbs->getLogP(AG1,GT2);break;
		}
	      break;
	    case 2:
	      switch(newPhase)
		{
		case 0:
		  transScore=transitionProbs->getLogP(AG2,GT0);break;
		case 1:
		  transScore=transitionProbs->getLogP(AG2,GT1);break;
		case 2:
		  transScore=transitionProbs->getLogP(AG2,GT2);break;
		}
	      break;
	    }
	  break;
	default: throw "bad! (1)";
	}
      break;
    case TAG:  
      switch(predType)
	{
	case ATG://                     * * * * * Single exon
	  transScore=transitionProbs->getLogP(ATG,TAG);break;
	case AG://                      * * * * * Final exon
	  switch(oldPhase)
	    {
	    case 0:
	      transScore=transitionProbs->getLogP(AG0,TAG);break;
	    case 1:
	      transScore=transitionProbs->getLogP(AG1,TAG);break;
	    case 2:
	      transScore=transitionProbs->getLogP(AG2,TAG);break;
	    }
	  break;
	default: throw "bad! (3)";
	}
      break;
    default: throw "bad! (2)";
    }
  return transScore;
}



void GeneZilla::observeRecentStopCodons(const BOOM::String &str,
				       SignalPtr signal)
{
  if(!prohibitPTCs) return;
  /*
    This method observes stop codons which occur anywhere from two bases 
    before the start of this signal's context window up to the last full 
    codon before the signal's consensus.  Any stop codons which are thus 
    found are stored in a set of phases which this signal cannot use for
    linking back to a coding signal.  This is necessary because at the
    time that this signal is instantiated, the main code for handling stop 
    codons will have processed only the sequence up to and including three
    bases left of the signal's context window, so any stop codons between
    that and the signal's consensus would otherwise be overlooked.
   */
  bool *p=recentlyEclipsedPhases;
  *p=*(p+1)=*(p+2)=false;
  int consPos=signal->getConsensusPosition();
  int begin=signal->getContextWindowPosition()-2;
  int end=consPos-3;
  switch(signal->getStrand()) {
  case FORWARD_STRAND:
    for(int pos=begin ; pos<=end ; ++pos)
      if(isochore->stopCodonSensor->consensusOccursAt(str,pos)) {
	int phase=posmod(consPos-pos);
	if(prohibitPTCs)
	  recentlyEclipsedPhases[phase]=true;
      }
    break;
  case REVERSE_STRAND:
    for(int pos=begin ; pos<=end ; ++pos)
      if(isochore->negStopCodonSensor->consensusOccursAt(str,pos)) {
	int phase=posmod(pos-consPos-1);
	if(prohibitPTCs)
	  recentlyEclipsedPhases[phase]=true;
      }
    break;
  }
}



int GeneZilla::mapPhaseBack(int phase,SignalPtr right,SignalPtr left)
{
  SignalType leftType=left->getSignalType();
  SignalType rightType=right->getSignalType();
  ContentType contentType=
    SignalTypeProperties::global.getContentType(leftType,rightType);
  switch(contentType)
    {
    case INITIAL_EXON:
    case INTERNAL_EXON:
    case FINAL_EXON:
    case SINGLE_EXON:
      {
	int length=
	  right->getConsensusPosition()-left->posOfBaseFollowingConsensus();
	return posmod(phase-length);
      }
    case INTRON:
    case UTR5_INTRON: 
    case UTR3_INTRON:
    case NEG_INTRON:
    case NEG_UTR5_INTRON:
    case NEG_UTR3_INTRON:
      return phase;
    case INTERGENIC:
      return (left->getStrand()==FORWARD_STRAND ? 0 : 2);
    case UTR5_INITIAL:
    case UTR5_INTERNAL:
    case UTR5_FINAL:
    case UTR5_SINGLE:
    case UTR3_INITIAL:
    case UTR3_INTERNAL:
    case UTR3_FINAL:
    case UTR3_SINGLE:
      return 0;
    case NEG_UTR5_INITIAL:
    case NEG_UTR5_INTERNAL:
    case NEG_UTR5_FINAL:
    case NEG_UTR5_SINGLE:
    case NEG_UTR3_INITIAL:
    case NEG_UTR3_INTERNAL:
    case NEG_UTR3_FINAL:
    case NEG_UTR3_SINGLE:
      return 2;
    case NEG_INITIAL_EXON:
    case NEG_INTERNAL_EXON:
    case NEG_FINAL_EXON:
    case NEG_SINGLE_EXON:
      {
	int length=
	  right->getConsensusPosition()-left->posOfBaseFollowingConsensus();
	return (phase+length)%3;
      }
    default: cout<<contentType<<endl; INTERNAL_ERROR;
    }
}



#ifdef EXPLICIT_GRAPHS
ParseGraph &GeneZilla::getParseGraph()
{
  return parseGraph;
}
#endif



DiscreteDistribution &GeneZilla::getDistribution(ContentType t)
{
  return *isochore->contentToDistribution[t];
}



ContentType GeneZilla::exonTypeToContentType(ExonType exonType,
					       Strand strand)
{
  return ::exonTypeToContentType(exonType,strand);
}



SignalSensor &GeneZilla::getSignalSensor(SignalType t)
{
  return *isochore->signalTypeToSensor[t];
}



ContentSensor &GeneZilla::getContentSensor(ContentType t)
{
  return *isochore->contentToSensor[t];
}



void GeneZilla::computeIntergenicSums(const Sequence &seq,
				     const BOOM::String &str,
				     const char *p)
{
  ContentSensor &sensor=getContentSensor(INTERGENIC);
  sensor.reset(seq,str,0);
  double score=0;
  for(int pos=0 ; pos<seqLen ; ++pos)
    {
      score+=sensor.scoreSingleBase(seq,str,pos,seq[pos],*p++);
      intergenicSums[pos]=score;
    }
  sensor.reset(seq,str,0);
}


#ifdef EXPLICIT_GRAPHS
BOOM::Vector<SignalPtr> *GeneZilla::getPathFromGff(
			       	  BOOM::Vector<BOOM::GffTranscript*> &t,
				  Sequence &seq,
				  const BOOM::String &seqStr,
				  BOOM::Vector<BOOL> &found)
{
  GffPathFromParseGraph g(*this);
  return g.getPathFromGff(t,seq,seqStr,found);
}
#endif


#ifdef FORCE_SPECIFIC_SIGNALS
void GeneZilla::forceSignalCoords(BOOM::Vector<int> &signalCoords)
{
  BOOM::Vector<int>::iterator cur=signalCoords.begin(), 
    end=signalCoords.end();
  for(; cur!=end ; ++cur)
    forcedSignalCoords.insert(*cur);
}
#endif



void GeneZilla::switchIsochore(float gcContent,int pos)
{
  // Change the active isochore
  isochore=isochores.getIsochore(gcContent);
  
  // Reset the content sensors (really only necessary for those implemented
  // as Finite State Machines)
  BOOM::Vector<ContentSensor*>::iterator cur=
    isochore->contentSensors.begin(),
    end=isochore->contentSensors.end();
  for(; cur!=end ; ++cur)
    {
      ContentSensor *sensor=*cur;
      sensor->reset(*seq,*seqStr,pos);
    }

  // Inform each signal queue of the change in isochore -- this just affects
  // which length distribution is used by the noncoding queues when they
  // are deciding which low-scoring signals to drop from their priority
  // queues -- and that only happens when EXPLICIT_GRAPHS are enabled
  BOOM::Vector<SignalQueue*>::iterator qCur=signalQueues.begin(),
    qEnd=signalQueues.end();
  for(; qCur!=qEnd ; ++qCur)
    {
      SignalQueue *queue=*qCur;
      queue->switchToIsochore(isochore);
    }
}



void GeneZilla::processIsochoreFile(const BOOM::String &filename,
				   float gcContent)
{
  isochores.load(filename);
  isochore=isochores.getIsochore(gcContent);
  createQueue(INITIAL_EXON);
  createQueue(INTERNAL_EXON);
  createQueue(FINAL_EXON);
  createQueue(SINGLE_EXON);
  createQueue(INTRON);
  createQueue(INTERGENIC);
  //createQueue(FIVE_PRIME_UTR);
  //createQueue(THREE_PRIME_UTR);
  createQueue(UTR5_INTRON);
  createQueue(UTR3_INTRON);
  createQueue(UTR5_INITIAL);
  createQueue(UTR5_INTERNAL);
  createQueue(UTR5_FINAL);
  createQueue(UTR5_SINGLE);
  createQueue(UTR3_INITIAL);
  createQueue(UTR3_INTERNAL);
  createQueue(UTR3_FINAL);
  createQueue(UTR3_SINGLE);

  switchIsochore(gcContent,0);
}



void GeneZilla::createQueue(ContentType contentType)
{
  ContentType revContentType=::reverseComplement(contentType);
  int queueCapacity=isochores.getQueueCapacity();
  SignalQueue *queue, *negQueue=NULL;
  NoncodingComparator *noncodingComp=
    isochore->noncodingComparators[contentType];
  BOOM::Array1D<SinglePhaseComparator*> *intronComp=
    isochore->intronComparators[contentType];
  switch(contentType)
    {
    case INTERGENIC:
      queue=new NoncodingQueue(contentType,queueCapacity,noncodingComp);
      signalQueues.push_back(queue);
      contentToQueue[contentType]=queue;
      break;
    case UTR5_INITIAL:
    case UTR5_INTERNAL:
    case UTR5_FINAL:
    case UTR5_SINGLE:
    case UTR3_INITIAL:
    case UTR3_INTERNAL:
    case UTR3_FINAL:
    case UTR3_SINGLE:
      queue=new NoncodingQueue(contentType,queueCapacity,noncodingComp);
      signalQueues.push_back(queue);
      contentToQueue[contentType]=queue;
#ifndef FORWARD_STRAND_ONLY
      negQueue=new NoncodingQueue(revContentType,queueCapacity,
				  noncodingComp);
      signalQueues.push_back(negQueue);
      contentToQueue[revContentType]=negQueue;
#endif
      break;

    case INTRON:
    case UTR5_INTRON:
    case UTR3_INTRON:
      queue=new IntronQueue(contentType,queueCapacity,*intronComp);
      signalQueues.push_back(queue);
      contentToQueue[contentType]=queue;
#ifndef FORWARD_STRAND_ONLY
      negQueue=new IntronQueue(revContentType,queueCapacity,*intronComp);
      signalQueues.push_back(negQueue);
      contentToQueue[revContentType]=negQueue;
#endif
      break;
    default:
      queue=new SignalQueue(contentType);
      signalQueues.push_back(queue);
      contentToQueue[contentType]=queue;
#ifndef FORWARD_STRAND_ONLY
      negQueue=new SignalQueue(revContentType);
      signalQueues.push_back(negQueue);
      contentToQueue[revContentType]=negQueue;
#endif
      break;
    }

  // Keep a separate list of coding queues (both strands separately)
  // so we can tell when a stop codon eclipses an ORF
  if(queue->precedesExon()) {
      forwardCodingQueues.push_back(queue);
#ifndef FORWARD_STRAND_ONLY
      reverseCodingQueues.push_back(negQueue);
#endif
    }

  // Register this queue with the appropriate content sensor in all 
  // isochores
  int numIsochores=isochores.getNumIsochores();
  for(int i=0 ; i<numIsochores ; ++i) {
    ContentSensor *sensor=
      isochores.getIthIsochore(i)->contentToSensor[contentType];
    sensor->addQueue(queue);

    if(negQueue) {
      sensor=
	isochores.getIthIsochore(i)->contentToSensor[revContentType];
      if(sensor) sensor->addQueue(negQueue);
      else cout<<revContentType<<" has a queue but no sensor!"<<endl;
    }
  }
}



void GeneZilla::resetAccumulatorPositions()
{
  BOOM::Vector<SignalQueue*>::iterator cur=signalQueues.begin(),
    end=signalQueues.end();
  for(; cur!=end ; ++cur)
    (*cur)->resetPosition();
}



void GeneZilla::loadIsochoreBoundaries(const BOOM::String &filename)
{
  IsochoreIntervalCmp cmp;
  BOOM::Map<BOOM::String,int> isoMap;
  isoMap["I"]=0;
  isoMap["II"]=1;
  isoMap["III"]=2;
  isoMap["IV"]=3;
  BOOM::GffReader reader(filename);
  BOOM::GffFeature *feature=NULL;
  while(feature=reader.nextFeature())
    {
      BOOM::String source=feature->getSource();
      source.tolower();
      if(source!="isochore") continue;
      const BOOM::String &isoSubstrate=feature->getSubstrate();
      int isochore=isoMap[feature->getFeatureType()];
      float GC=-1;
      BOOM::Vector<BOOM::String> &fields=feature->getExtraFields();
      int n=fields.size();
      for(int i=0 ; i<n ; ++i)
	{
	  BOOM::String field=fields[i];
	  if(GCregex.search(field))
	    GC=GCregex[1].asInt()/100.0;
	}
      if(GC<0) throw 
	BOOM::String("Could not parse \"gc\" field from GFF record:\n")+
	feature->getRawLine();
      isochoreIntervals[isoSubstrate].push_back(
	    IsochoreInterval(feature->getBegin(),
			     feature->getEnd(),
			     isochore,GC));
      delete feature;
    }

  // Sort boundaries for each substrate
  BOOM::Map< BOOM::String,BOOM::Vector<IsochoreInterval> >::iterator
    cur=isochoreIntervals.begin(), end=isochoreIntervals.end();
  for(; cur!=end ; ++cur)
    {
      BOOM::Vector<IsochoreInterval> &boundaries=(*cur).second;
      VectorSorter<IsochoreInterval> sorter(boundaries,cmp);
      sorter.sortAscendInPlace();
    }
}



void GeneZilla::crossIsochoreBoundary(int pos)
{
  gcContent=nextIsochoreInterval.GC;
  switchIsochore(gcContent,pos);
  ++nextIsochoreIndex;
  BOOM::Vector<IsochoreInterval> &intervals=isochoreIntervals[substrateId];
  if(nextIsochoreIndex<intervals.size())
    nextIsochoreInterval=intervals[nextIsochoreIndex];
  else
    nextIsochoreInterval.begin=-1;
  cerr<<pos<<" switching to isochore "<<isochore->name<<endl;
}



void GeneZilla::loadCpGislands(const BOOM::String &filename)
{
  throw "CpG state not yet implemented";
}



BOOM::StringMap<char> *GeneZilla::getStopCodonConsensuses()
{
  return 
    &isochores.getIthIsochore(0)->stopCodonSensor->getConsensuses();
}



Transitions *GeneZilla::getTransitionProbs()
{
  return isochore->transitionProbs;
}



double GeneZilla::scoreExon(SignalPtr thisSignal,SignalPtr nextSignal,
			   int phase,ContentType &contentType)
{
  SignalType thisType=thisSignal->getSignalType();
  SignalType nextType=nextSignal->getSignalType();
  int begin=thisSignal->posOfBaseFollowingConsensus();
  int end=nextSignal->getConsensusPosition();
  int length=end-begin;
  int nextPhase=getPhase(nextSignal);
  contentType=
    SignalTypeProperties::global.getContentType(thisType,nextType);
  double lengthPenalty=useDurationScores ? 
    isochore->contentToDistribution[contentType]->getLogP(length) : 0;
  double numerator=
    nextSignal->getInductiveScore(nextPhase)-
    nextSignal->contextWindowScore() -
    thisSignal->getInductiveScore(phase) -
    isochore->transitionProbs->getLogP(thisType,nextType) -
    lengthPenalty;
  double denominator=
    intergenicSums[begin+length-1]-(begin ? intergenicSums[begin-1] : 0);
  return numerator-denominator;
}



void GeneZilla::generateGff(BOOM::Stack<SignalPtr> *path,int seqLen,
			   double parseScore)
{
  // Invariant: the stack always contains at least 2 signals, and both
  // the top signal and the bottom signal are bogus (left & right termini)
  SignalPtr thisSignal=path->pop();
  SignalPtr firstSignal=thisSignal; // should be LEFT_TERMINUS

  setPrecision(cout,2);
  cout << "Single best splice pattern:" << endl;

  while(!path->isEmpty()) {
    int phase=getPhase(thisSignal);
    SignalType thisType=thisSignal->getSignalType();
    Strand strand=thisSignal->getStrand();
    SignalPtr nextSignal=path->pop();
    SignalType nextType=nextSignal->getSignalType();
    int thisPos=thisSignal->getConsensusPosition();
    int nextPos=nextSignal->getConsensusPosition();
    if(thisPos<0) thisPos=0;
    double exonScore;

    // Compute exon score
    ContentType contentType;
    if(::beginsCoding(thisType)) 
      exonScore=scoreExon(thisSignal,nextSignal,phase,contentType);

    switch(thisType){
    case TSS:
      cout << substrateId << '\t' << PROGRAM_NAME << "\texon\t"
	   << thisPos+1 << '\t' 
	   << nextPos
	   << "\t.\t"<<strand<<"\t.\t"
	   << "transcript_id=" 
	   << (thisType==TSS ? transcriptId+1 : transcriptId)
	   << ";" << endl;
      break;
    case UTR5AG:
      cout << substrateId<<'\t'<<PROGRAM_NAME<<"\texon\t"
	   << thisPos+3 << '\t' 
	   << nextPos
	   << "\t.\t"<<strand<<"\t.\t"
	   << "transcript_id=" 
	   << (thisType==TES ? transcriptId : transcriptId+1)
	   << ";" << endl;
      break;
    }
    thisSignal=nextSignal;
  }
  
  cout<<"Sequence length: " << seqLen << endl;
  cout<<"Parse score: " << parseScore << endl;
}



