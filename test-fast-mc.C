/****************************************************************
 test-fast-mc
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <math.h>
#include <iostream>
#include <fstream>
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/ConfigFile.H"
#include "BOOM/Map.H"
#include "BOOM/Stack.H"
#include "BOOM/Constants.H"
#include "BOOM/Time.H"
#include "genezilla.H"
#include "SignalSensor.H"
#include "SignalQueue.H"
#include "SignalTypeProperties.H"
#include "NoncodingQueue.H"
#include "IntronQueue.H"
#include "EmpiricalDistribution.H"
#include "GeometricDistribution.H"
#include "Transitions.H"
#include "FastMarkovChain.H"

static const char *PROGRAM_NAME="test-fast-mc";
static const char *VERSION="1.0";
Alphabet alphabet;
int frame; // ### CAUTION: this is required by older code; to be removed

class GeneZilla
{
  BOOM::String substrateId;
  BOOM::Vector<SignalSensor*> signalSensors;
  BOOM::Vector<SignalQueue*> signalQueues;
  BOOM::Vector<SignalQueue*> forwardCodingQueues, reverseCodingQueues;
  BOOM::Map<ContentType,SignalQueue*> contentToQueue;
  SignalSensor *stopCodonSensor, *negStopCodonSensor;
  Transitions *transitionProbs;

  float getGCcontent(BOOM::String &);
  void processConfigFile(BOOM::ConfigFile &);
  void loadSubmodels(BOOM::ConfigFile &);
  void loadSignalSensor(BOOM::ConfigFile &,const BOOM::String &modelLabel,
			const BOOM::String &consensusLabel);
  void loadContentSensor(BOOM::ConfigFile &,const BOOM::String &modelLabel,
			 DistributionType,const BOOM::String &distrLabel);
  void mainAlgorithm(Sequence &,const BOOM::String &);
  void instantiateLeftTermini();
  BOOM::Stack<Signal*> *instantiateRightTermini(int seqLen);
  inline void updateAccumulators(Sequence &,const BOOM::String &,int pos,
				 Symbol,char);
  /*inline*/ void linkBack(Signal &newSignal);
  inline void selectPredecessors(int newConsPos,SignalQueue &queue,
				 ContentType contentType,Strand strand,
				 double bestScore[3],Signal *bestPred[3],
				 SignalType toType);
  inline void selectIntergenicPred(int newConsPos,SignalQueue &queue,
				   Strand strand,double bestScore[3],
				   Signal *bestPred[3],ContentType,
				   SignalType toType);
  inline void selectCodingPred(int newConsPos,SignalQueue &queue,Strand strand,
			       double bestScore[3],Signal *bestPred[3],
			       ContentType,SignalType toType);
  inline void selectIntronPred(int newConsPos,SignalQueue &queue,Strand strand,
			       double bestScore[3],Signal *bestPred[3],
			       ContentType,SignalType toType);
  void enqueue(Signal &);
  void handleStopCodons(const BOOM::String &,int pos);
  void terminateForwardORFs(int position);
  void terminateReverseORFs(int position);
  void loadTransProbs(const BOOM::String &transFile);
  BOOM::Stack<Signal*> *traceBack(Signal *rightTerminus,int phase);
  void generateGff(BOOM::Stack<Signal*> *path);
public:
  GeneZilla();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	GeneZilla app;
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



GeneZilla::GeneZilla()
  {
    // ctor
  }



int GeneZilla::main(int argc,char *argv[])
  {
    // Process command line
    BOOM::CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=2)
      throw string("test-fast-mc <*.cfg> <*.fasta>");
    BOOM::String configFilename=cmd.arg(0);
    BOOM::String fastaFilename=cmd.arg(1);

    cout << "Loading sequence..." << endl;
    alphabet=DnaAlphabet::global;
    Sequence *substrate=Sequence::load(fastaFilename,alphabet,substrateId);
    BOOM::String *substrateStr=substrate->toString(alphabet);
    Sequence *revSubstrate=substrate->reverseComplement(alphabet);
    BOOM::String *revSubstrateStr=revSubstrate->toString(alphabet);

    cout << "Processing config file..." << endl;
    BOOM::ConfigFile configFile(configFilename);
    processConfigFile(configFile);

    ContentSensor &sensor=
      contentToQueue[INTERNAL_EXON]->getContentSensor();
    ContentSensor &revSensor=*sensor.reverseComplement();

    int len=substrate->getLength();
    for(int j=0 ; j<12 ; ++j)
      {
	/* TEST 1: reverse-complementation works for WMM!
	double fwScore=
	  stopCodonSensor->getLogP(*substrate,*substrateStr,j);
	double revScore=
	  negStopCodonSensor->getLogP(*revSubstrate,*revSubstrateStr,
				      len-j-stopCodonSensor->
				        getContextWindowLength());
	*****************************************************************/

	/* TEST 2: reverse-complementation works for simple MarkovChain!
	double fwScore=
	  sensor.scoreSingleBase(*substrate,*substrateStr,j,(*substrate)[j],
				 (*substrateStr)[j]);
	double revScore=
	  revSensor.scoreSingleBase(*revSubstrate,*revSubstrateStr,len-1-j,
				 (*revSubstrate)[len-1-j],
				 (*revSubstrateStr)[len-1-j]);
	******************************************************************/


	/* TEST 3: reverse-complementation works for ThreePeriodicMarkovChain!

	double fwPhase0, fwPhase1, fwPhase2;
	sensor.scoreSingleBase(*substrate,*substrateStr,j,(*substrate)[j],
			       (*substrateStr)[j],fwPhase0,fwPhase1,fwPhase2);
	double revPhase0, revPhase1, revPhase2;
	revSensor.scoreSingleBase(*substrate,*substrateStr,j,
				  (*substrate)[j],
				  (*substrateStr)[j],
				  revPhase0,revPhase1,revPhase2);

	cout << j << " (0) " << fwPhase0 << " vs. " << revPhase0 << endl;
	cout << "   (1) " << fwPhase1 << " vs. " << revPhase1 << endl;
	cout << "   (2) " << fwPhase2 << " vs. " << revPhase2 << endl;
	******************************************************************/

      }

    return 0;
  }



float GeneZilla::getGCcontent(BOOM::String &seq)
{
  int n=seq.length(), ATCG=0, GC=0;
  const char *p=seq.c_str();
  for(int i=0 ; i<n ; ++i)
    switch(*p++)
      {
      case 'G':
      case 'C':
	++GC;
      case 'A':
      case 'T':
	++ATCG;
      }
  return GC/float(ATCG);
}



void GeneZilla::loadSubmodels(BOOM::ConfigFile &configFile)
{
  loadSignalSensor(configFile,"donor-model","donor-consensus");
  loadSignalSensor(configFile,"acceptor-model","acceptor-consensus");
  loadSignalSensor(configFile,"start-codon-model","start-codon-consensus");
  loadSignalSensor(configFile,"stop-codon-model","stop-codon-consensus");
  loadSignalSensor(configFile,"polya-model","polya-consensus");
  loadSignalSensor(configFile,"promoter-model","promoter-consensus");

  loadContentSensor(configFile,"initial-exons",EMPIRICAL_DISTRIBUTION,
		    "initial-exon-lengths");
  loadContentSensor(configFile,"internal-exons",EMPIRICAL_DISTRIBUTION,
		    "internal-exon-lengths");
  loadContentSensor(configFile,"final-exons",EMPIRICAL_DISTRIBUTION,
		    "final-exon-lengths");
  loadContentSensor(configFile,"single-exons",EMPIRICAL_DISTRIBUTION,
		    "single-exon-lengths");
  loadContentSensor(configFile,"introns",GEOMETRIC_DISTRIBUTION,
		    "mean-intron-length");
  loadContentSensor(configFile,"intergenic",GEOMETRIC_DISTRIBUTION,
		    "mean-intergenic-length");
  loadContentSensor(configFile,"3'-UTR",GEOMETRIC_DISTRIBUTION,
		    "mean-3'-UTR-length");
  loadContentSensor(configFile,"5'-UTR",GEOMETRIC_DISTRIBUTION,
		    "mean-5'-UTR-length");
}



void GeneZilla::loadContentSensor(BOOM::ConfigFile &configFile,
				 const BOOM::String &modelLabel,
				 DistributionType distrType,
				 const BOOM::String &distrLabel)
{
  // Load the model
  BOOM::String filename=configFile.lookupOrDie(modelLabel);
  ContentSensor *sensor=ContentSensor::load(filename);
  
  // Load or create the distribution
  DiscreteDistribution *distribution;
  BOOM::String distrParm=configFile.lookupOrDie(distrLabel);
  if(distrType==EMPIRICAL_DISTRIBUTION)
    distribution=new EmpiricalDistribution(distrParm);
  else
    distribution=new GeometricDistribution(distrParm.asInt());

  // Create a signal queue for it
  ContentType contentType=sensor->getContentType();
  SignalQueue *queue, *negQueue=NULL;
  switch(contentType)
    {
    case INTERGENIC:
      queue=new NoncodingQueue(*sensor,*distribution);
      signalQueues.push_back(queue);
      contentToQueue[contentType]=queue;
      break;
    case FIVE_PRIME_UTR:
    case THREE_PRIME_UTR:
    //case NEG_FIVE_PRIME_UTR:
    //case NEG_THREE_PRIME_UTR:
      queue=new NoncodingQueue(*sensor,*distribution);
      negQueue=new NoncodingQueue(*sensor->reverseComplement(),
				  *distribution);
      signalQueues.push_back(queue);
      signalQueues.push_back(negQueue);
      contentToQueue[contentType]=queue;
      contentToQueue[reverseComplement(contentType)]=negQueue;
      break;
    case INTRON:
    //case NEG_INTRON:
      queue=new IntronQueue(*sensor,*distribution);
      negQueue=new IntronQueue(*sensor->reverseComplement(),*distribution);
      signalQueues.push_back(queue);
      signalQueues.push_back(negQueue);
      contentToQueue[contentType]=queue;
      contentToQueue[reverseComplement(contentType)]=negQueue;
      break;
    default:
      queue=new SignalQueue(*sensor,*distribution);
      negQueue=new SignalQueue(*sensor->reverseComplement(),*distribution);
      signalQueues.push_back(queue);
      signalQueues.push_back(negQueue);
      contentToQueue[contentType]=queue;
      contentToQueue[reverseComplement(contentType)]=negQueue;
      break;
    }

  // Keep a separate list of coding queues (both strands separately)
  // so we can tell when a stop codon eclipses an ORF
  if(queue->precedesExon())
    {
      forwardCodingQueues.push_back(queue);
      reverseCodingQueues.push_back(negQueue); // ### is this right?
    }
}



void GeneZilla::loadSignalSensor(BOOM::ConfigFile &configFile,
				   const BOOM::String &modelLabel,
				   const BOOM::String &consensusLabel)
{
  // Load the model
  BOOM::String filename=configFile.lookupOrDie(modelLabel);
  SignalSensor *sensor=SignalSensor::load(filename);
  
  // Load consensus strings
  BOOM::String consensuses=configFile.lookupOrDie(consensusLabel);
  BOOM::Vector<BOOM::String> *fields=consensuses.getFields(",|;:/.&_-+");
  BOOM::Vector<BOOM::String>::iterator cur=fields->begin(), end=fields->end();
  for(; cur!=end ; ++cur)
    sensor->addConsensus(*cur);
  delete fields;

  // Reverse-complement the model
  SignalSensor *negSensor=sensor->reverseComplement();

  // Add model to global model list
  signalSensors.push_back(sensor);
  signalSensors.push_back(negSensor);

  if(sensor->getSignalType()==TAG)
    {
      stopCodonSensor=sensor;
      negStopCodonSensor=negSensor;
    }
}



void GeneZilla::processConfigFile(BOOM::ConfigFile &configFile)
{
  // The following block will eventually be moved into the 
  // configuration file, which will allow user-defined state
  // topologies:
  SignalTypeProperties &stp=SignalTypeProperties::global;
  stp.allowPhase(ATG,0);
  stp.allowPhase(TAG,0);
  stp.allowPhase(PROM,0);
  stp.allowPhase(POLYA,0);
  stp.allowPhase(GT,0).allowPhase(GT,1).allowPhase(GT,2);
  stp.allowPhase(AG,0).allowPhase(AG,1).allowPhase(AG,2);
  stp.allowPhase(NEG_ATG,2);
  stp.allowPhase(NEG_TAG,2);
  stp.allowPhase(NEG_PROM,2);
  stp.allowPhase(NEG_POLYA,2);
  stp.allowPhase(NEG_GT,0).allowPhase(NEG_GT,1).allowPhase(NEG_GT,2);
  stp.allowPhase(NEG_AG,0).allowPhase(NEG_AG,1).allowPhase(NEG_AG,2);
  stp.belongTo(ATG,INITIAL_EXON).belongTo(ATG,SINGLE_EXON);
  stp.belongTo(TAG,THREE_PRIME_UTR).belongTo(TAG,INTERGENIC);
  stp.belongTo(GT,INTRON);
  stp.belongTo(AG,INTERNAL_EXON).belongTo(AG,FINAL_EXON);
  stp.belongTo(PROM,FIVE_PRIME_UTR);
  stp.belongTo(POLYA,INTERGENIC);
  stp.belongTo(NEG_ATG,INTERGENIC).belongTo(NEG_ATG,NEG_FIVE_PRIME_UTR);
  stp.belongTo(NEG_TAG,NEG_SINGLE_EXON).belongTo(NEG_TAG,NEG_FINAL_EXON);
  stp.belongTo(NEG_GT,NEG_INITIAL_EXON).belongTo(NEG_GT,NEG_INTERNAL_EXON);
  stp.belongTo(NEG_AG,NEG_INTRON);
  stp.belongTo(NEG_PROM,INTERGENIC);
  stp.belongTo(NEG_POLYA,NEG_THREE_PRIME_UTR);
  stp.linkTo(NEG_ATG,NEG_INITIAL_EXON).linkTo(NEG_ATG,NEG_SINGLE_EXON);
  stp.linkTo(NEG_TAG,NEG_THREE_PRIME_UTR).linkTo(NEG_TAG,INTERGENIC);
  stp.linkTo(NEG_GT,NEG_INTRON);
  stp.linkTo(NEG_AG,NEG_INTERNAL_EXON).linkTo(NEG_AG,NEG_FINAL_EXON);
  stp.linkTo(NEG_PROM,NEG_FIVE_PRIME_UTR);
  stp.linkTo(NEG_POLYA,INTERGENIC);
  stp.linkTo(ATG,INTERGENIC).linkTo(ATG,FIVE_PRIME_UTR);
  stp.linkTo(TAG,SINGLE_EXON).linkTo(TAG,FINAL_EXON);
  stp.linkTo(GT,INITIAL_EXON).linkTo(GT,INTERNAL_EXON);
  stp.linkTo(AG,INTRON);
  stp.linkTo(PROM,INTERGENIC);
  stp.linkTo(POLYA,THREE_PRIME_UTR);
  stp.setStrand(ATG,FORWARD_STRAND);
  stp.setStrand(GT,FORWARD_STRAND);
  stp.setStrand(AG,FORWARD_STRAND);
  stp.setStrand(TAG,FORWARD_STRAND);
  stp.setStrand(PROM,FORWARD_STRAND);
  stp.setStrand(POLYA,FORWARD_STRAND);
  stp.setStrand(NEG_ATG,REVERSE_STRAND);
  stp.setStrand(NEG_GT,REVERSE_STRAND);
  stp.setStrand(NEG_AG,REVERSE_STRAND);
  stp.setStrand(NEG_TAG,REVERSE_STRAND);
  stp.setStrand(NEG_PROM,REVERSE_STRAND);
  stp.setStrand(NEG_POLYA,REVERSE_STRAND);
  stp.setConsensusCoding(ATG,true);
  stp.setConsensusCoding(GT,false);
  stp.setConsensusCoding(AG,false);
  stp.setConsensusCoding(TAG,true);
  stp.setConsensusCoding(PROM,false);
  stp.setConsensusCoding(POLYA,false);
  stp.setConsensusCoding(NEG_ATG,true);
  stp.setConsensusCoding(NEG_GT,false);
  stp.setConsensusCoding(NEG_AG,false);
  stp.setConsensusCoding(NEG_TAG,true);
  stp.setConsensusCoding(NEG_PROM,false);
  stp.setConsensusCoding(NEG_POLYA,false);

  cout << "Loading submodels..." << endl;
  loadSubmodels(configFile);

  cout << "Loading transition probabilities..." << endl;
  BOOM::String transFile=configFile.lookupOrDie("transition-probabilities");
  loadTransProbs(transFile);
}



void GeneZilla::loadTransProbs(const BOOM::String &transFile)
{
  ifstream is(transFile.c_str());
  if(!is.good()) throw BOOM::String("Can't open file ")+transFile;
  int numSignalTypes;
  is >> numSignalTypes;
  transitionProbs=new Transitions(numSignalTypes,is);
}



void GeneZilla::instantiateLeftTermini()
{
  BOOM::Vector<SignalSensor*>::iterator cur=signalSensors.begin(),
    end=signalSensors.end();
  for(; cur!=end ; ++cur)
    {
      SignalSensor &s=**cur;
      BOOM::Set<ContentType> &queues=s.belongsInWhichQueues();
      BOOM::Set<ContentType>::iterator cur=queues.begin(), end=queues.end();
      for(; cur!=end ; ++cur)
	{
	  SignalQueue *queue=contentToQueue[*cur];
#ifdef DEBUG
	  cout<<"considering adding "<<s.getSignalType()<<" to "<<queue->getContentType()<<" queue..."<<endl;
#endif
	  if(queue->isEmpty())
	    {
	      queue->addSignal(s.getLeftTerminus());
#ifdef DEBUG
	      cout<<"  ===>adding "<<s.getSignalType()<<" to "<<queue->getContentType()<<" queue..."<<endl;
#endif
	    }
	}
    }
}



BOOM::Stack<Signal*> *GeneZilla::instantiateRightTermini(int seqLen)
{
  int bestPhase;
  double bestScore;
  Signal *bestSignal;
  bestScore=NEGATIVE_INFINITY;
  bestSignal=NULL;

  BOOM::Set<ContentType> instantiated;
  BOOM::Vector<SignalSensor*>::iterator cur=signalSensors.begin(),
    end=signalSensors.end();
  for(; cur!=end ; ++cur)
    {
      SignalSensor &s=**cur;
      BOOM::Set<ContentType> &queues=s.belongsInWhichQueues();
      BOOM::Set<ContentType>::iterator cur=queues.begin(), end=queues.end();
      for(; cur!=end ; ++cur)
	{
	  ContentType contentType=*cur;
	  if(!instantiated.isMember(contentType))
	    {
	      SignalQueue *queue=contentToQueue[contentType];
	      Signal *rightTerminus=s.getRightTerminus(seqLen);
	      queue->addSignal(rightTerminus);
	      linkBack(*rightTerminus);
	      Propagator &prop=rightTerminus->getPropagator(contentType);
	      for(int i=0 ; i<3 ; ++i)
		{
		  double thisScore=prop[i];
		  if(thisScore>bestScore)
		    {
		      bestScore=thisScore;
		      bestSignal=rightTerminus;
		      bestPhase=i;
		    }
		}
	    }
	}
    }

  return traceBack(bestSignal,bestPhase);
}



void GeneZilla::mainAlgorithm(Sequence &seq,const BOOM::String &str)
{
  // Instantiate one signal of each type at the left terminus to act as
  // anchors to which real signals can link back
  instantiateLeftTermini();

  // Make a single left-to-right pass across the sequence
  int seqLen=seq.getLength();
  const char *charPtr=str.c_str();
  for(int pos=0 ; pos<seqLen ; ++pos, ++charPtr)
    {
      //cout<<"getting next base"<<endl;
      Symbol base=seq[pos];
      
      // Check whether any signals occur here
      BOOM::Vector<SignalSensor*>::iterator cur=signalSensors.begin(),
	end=signalSensors.end();
      for(; cur!=end ; ++cur )
	{
	  SignalSensor &sensor=**cur;
	  //cout << "sensor.detect()" << endl;
	  Signal *signal=sensor.detect(seq,str,pos);
	  if(signal)
	    {
	      // Find optimal predecessor for this signal in all 3 phases
	      linkBack(*signal);

	      // Add this signal to the appropriate queue(s)
	      enqueue(*signal);

	      // If this is a stop codon, terminate this reading frame
	      handleStopCodons(str,pos);
	    }
	}
      
      // Propagate scores of all non-eclipsed signals up to this base
      updateAccumulators(seq,str,pos,base,*charPtr);
    }

  // Instantiate an anchor signal of each type at the right terminus
  // and link them back, to complete the dynamic programming evaluation:
  BOOM::Stack<Signal*> *path=instantiateRightTermini(seqLen);
  generateGff(path);
}



inline void GeneZilla::updateAccumulators(Sequence &seq,const BOOM::String &str,
					 int pos,Symbol base,char c)
{
  BOOM::Vector<SignalQueue*>::iterator cur=signalQueues.begin(),
    end=signalQueues.end();
  for(; cur!=end ; ++cur)
    {
      SignalQueue &queue=**cur;
      ContentSensor &contentSensor=queue.getContentSensor();
      if(contentSensor.isCoding())
	{
	  double scorePhase0, scorePhase1, scorePhase2;
	  contentSensor.scoreSingleBase(seq,str,pos,base,c,scorePhase0,
					scorePhase1,scorePhase2);
	  queue.addToAccumulator(scorePhase0,scorePhase1,scorePhase2,pos);
	  //if(queue.getContentType()==FIVE_PRIME_UTR) cout<<pos<<":"<<queue.getContentType()<<"=("<<scorePhase0<<","<<scorePhase1<<","<<scorePhase2<<") : "<<queue.getAccumulator()<<endl;
	}
      else
	{
	  double score=contentSensor.scoreSingleBase(seq,str,pos,base,c);
	  queue.addToAccumulator(score);
	  //if(queue.getContentType()==INTRON) cout<<pos<<":"<<queue.getContentType()<<"="<<score<<" : "<<queue.getAccumulator()<<endl;
	}
    }
}



/*inline*/ void GeneZilla::linkBack(Signal &newSignal)
{
#ifdef DEBUG
  cout << "LINKING BACK " << newSignal.getSignalType() << " @ " << newSignal.getConsensusPosition() <<" "<<newSignal.contextWindowScore()<< endl;
#endif
  Signal *bestPred[3];
  bestPred[0]=bestPred[1]=bestPred[2]=NULL;
  double bestScore[3];
  bestScore[0]=bestScore[1]=bestScore[2]=NEGATIVE_INFINITY;
  int newConsPos=newSignal.getConsensusPosition();
  Strand strand=newSignal.getStrand();
  SignalType signalType=newSignal.getSignalType();

  // Consider all queues that newSignal could link back to
  BOOM::Set<ContentType> &queues=newSignal.linksBackToWhichQueues();
  BOOM::Set<ContentType>::iterator cur=queues.begin(), end=queues.end();
  for(; cur!=end ; ++cur)
    {
      // Make sure the queue's signals are all propagated up to this point
      ContentType contentType=*cur;
#ifdef DEBUG
      cout << "contentType=" << contentType << endl;
#endif
      if(!contentToQueue.isDefined(contentType)) throw "bad!";//###DEBUGGING
      SignalQueue *queue=contentToQueue[contentType];
      queue->flushAccumulator();
#ifdef DEBUG
      cout << "\t\t|Q|="<<queue->numElements()<<endl;
#endif

      // Consider all of the queue's signals as potential predecessors
      selectPredecessors(newConsPos,*queue,contentType,strand,
			 bestScore,bestPred,signalType);
    }
  // Install the selected predecessors as dyn. prog. links and update
  // inductive scores for all propagators of this signal (identically)
  BOOM::Set<ContentType> &nextQueues=newSignal.belongsInWhichQueues();
  cur=nextQueues.begin();
  end=nextQueues.end();
  for(; cur!=end ; ++cur)
    {
      ContentType contentType=*cur;
      Propagator &prop=newSignal.getPropagator(contentType);
#ifdef DEBUG
      cout<<"installing into slot "<<(SignalTypeProperties::global.whichPropagator(newSignal.getSignalType(),contentType))<<" for "<<newSignal.getSignalType()<<" x "<<contentType<<endl;
      cout<<"before installation, prop="<<prop<<endl;
#endif
      prop[0]+=bestScore[0];
      prop[1]+=bestScore[1];
      prop[2]+=bestScore[2];
#ifdef DEBUG
      cout<<"installing "<<bestScore[0]<<","<<bestScore[1]<<","<<bestScore[2]<<" for "<<contentType<<"=>"<<prop<<endl;
#endif
    }
  newSignal.setPredecessor(0,bestPred[0]);
  newSignal.setPredecessor(1,bestPred[1]);
  newSignal.setPredecessor(2,bestPred[2]);
}



inline void GeneZilla::selectPredecessors(int newConsPos,SignalQueue &queue,
					 ContentType contentType,Strand strand,
					 double bestScore[3],
					 Signal *bestPred[3],
					 SignalType toType)
{
  switch(contentType)
    {
    case FIVE_PRIME_UTR:
    case THREE_PRIME_UTR:
    case NEG_FIVE_PRIME_UTR:
    case NEG_THREE_PRIME_UTR:
    case INTERGENIC:
      selectIntergenicPred(newConsPos,queue,strand,bestScore,bestPred,
			   contentType,toType);
      break;
    case INTRON:
    case NEG_INTRON:
      selectIntronPred(newConsPos,queue,strand,bestScore,bestPred,
		       contentType,toType);
      break;
    default:
      selectCodingPred(newConsPos,queue,strand,bestScore,bestPred,
		       contentType,toType);
    }
}



inline void GeneZilla::selectCodingPred(int newConsPos,SignalQueue &queue,
				       Strand strand,double bestScore[3],
				       Signal *bestPred[3],
				       ContentType contentType,
				       SignalType toType)
{
  // Version #1 of 3: CODING.  Phase is translated across the exon.

  SignalQueueIterator cur=queue.begin(), end=queue.end();
  for(; cur!=end ; ++cur)
    {
      Signal &pred=**cur;
      int oldPos=pred.posOfBaseFollowingConsensus();
      int length=newConsPos-oldPos;
      double lengthScore=queue.getDistribution().getLogP(length);
      double transScore=transitionProbs->getLogP(pred.getSignalType(),toType);
#ifdef DEBUG
      cout<<"length="<<length<<" ("<<oldPos<<" to "<<newConsPos<<") "<<" lenscore="<<lengthScore<<" transscr="<<transScore<<"("<<pred.getSignalType()<<"->"<<toType<<")"<<endl;
#endif
      int frameDelta=length % 3;
#ifdef DEBUG
      cout<<"frameDelta="<<frameDelta<<endl;
#endif
      Propagator &predProp=pred.getPropagator(contentType);
      for(int oldPhase=0 ; oldPhase<3 ; ++oldPhase)
	{
	  int newPhase=(strand==FORWARD_STRAND ?
			(oldPhase+frameDelta) % 3 :
			posmod(oldPhase-frameDelta));
	  double predScore=predProp[oldPhase] + lengthScore + transScore;
	  double &bestPredScore=bestScore[newPhase];
	  if(finite(predScore) && predScore>bestPredScore)
	    {
	      bestPredScore=predScore;
	      bestPred[newPhase]=&pred;
#ifdef DEBUG
	      cout << "\t\tselected predecessor for phase " << newPhase<<", score="<<predScore<<endl;
#endif
	    }
	}
    }
}



inline void GeneZilla::selectIntronPred(int newConsPos,SignalQueue &queue,
				       Strand strand,double bestScore[3],
				       Signal *bestPred[3],
				       ContentType contentType,
				       SignalType toType)
{
  // Version #2 of 3: INTRONS.  Phase is not translated across the intron.

  SignalQueueIterator cur=queue.begin(), end=queue.end();
  for(; cur!=end ; ++cur)
    {
      Signal &pred=**cur;
      double transScore=transitionProbs->getLogP(pred.getSignalType(),toType);
      int oldPos=pred.posOfBaseFollowingConsensus();
      int length=newConsPos-oldPos;
      double lengthScore=queue.getDistribution().getLogP(length);
#ifdef DEBUG
      cout<<"length="<<length<<" ("<<oldPos<<" to "<<newConsPos<<") "<<" lenscore="<<lengthScore<<" transscr="<<transScore<<"("<<pred.getSignalType()<<"->"<<toType<<")"<<endl;
#endif
      Propagator &predProp=pred.getPropagator(contentType);
      for(int phase=0 ; phase<3 ; ++phase)
	{
	  double predScore=predProp[phase] + lengthScore + transScore;
	  double &bestPredScore=bestScore[phase];
	  if(finite(predScore) && predScore>bestPredScore)
	    {
	      bestPredScore=predScore;
	      bestPred[phase]=&pred;
#ifdef DEBUG
	      cout << "\t\tselected predecessor for phase " << phase<< endl;
#endif
	    }
	}
    }
}



inline void GeneZilla::selectIntergenicPred(int newConsPos,SignalQueue &queue,
					   Strand strand,double bestScore[3],
					   Signal *bestPred[3],
					   ContentType contentType,
					   SignalType toType)
{
  // Version #3 of 3: INTERGENIC.  Phase 0 or 2 is used depending on strand.

  int newPhase=(strand==FORWARD_STRAND ? 0 : 2);
  SignalQueueIterator cur=queue.begin(), end=queue.end();
  for(; cur!=end ; ++cur)
    {
      Signal &pred=**cur;
#ifdef DEBUG
      cout<<"\tpotential pred="<<pred.getSignalType()<<"@"<<pred.getContextWindowPosition()<<endl;
#endif
      double transScore=transitionProbs->getLogP(pred.getSignalType(),toType);
      int oldPos=pred.posOfBaseFollowingConsensus();
      int length=newConsPos-oldPos;
      double lengthScore=queue.getDistribution().getLogP(length);
#ifdef DEBUG
      cout<<"length="<<length<<" ("<<oldPos<<" to "<<newConsPos<<") "<<" lenscore="<<lengthScore<<" transscr="<<transScore<<" "<<queue.getContentType()<<endl;
#endif
      Propagator &predProp=pred.getPropagator(contentType);
      int oldPhase=(pred.getStrand()==FORWARD_STRAND ? 0 : 2);
      double predScore=predProp[oldPhase] + lengthScore + transScore;
      double &bestPredScore=bestScore[newPhase];
      //cout << "\t\t\t" << pred.getSignalType()<<" "<<predScore<<" "<<oldPhase<<" "<<bestPredScore<<" "<<newPhase<<endl;
      if(finite(predScore) && predScore>bestPredScore)
	{
	  bestPredScore=predScore;
	  bestPred[newPhase]=&pred;
#ifdef DEBUG
	  cout << "\t\tselected predecessor for phase " << newPhase<< endl;
#endif
	}
    }
}



void GeneZilla::enqueue(Signal &signal)
{
  BOOM::Set<ContentType> &queues=signal.belongsInWhichQueues();
  BOOM::Set<ContentType>::iterator cur=queues.begin(), end=queues.end();
  for(; cur!=end ; ++cur)
    {
      ContentType contentType=*cur;
      contentToQueue[contentType]->addSignal(&signal);
    }
}



void GeneZilla::handleStopCodons(const BOOM::String &str,int pos)
{
  if(stopCodonSensor->consensusOccursAt(str,pos))
    terminateForwardORFs(pos);
  else if(negStopCodonSensor->consensusOccursAt(str,pos))
    terminateReverseORFs(pos);
}



void GeneZilla::terminateForwardORFs(int TAGposition)
{
  BOOM::Vector<SignalQueue*>::iterator qCur=forwardCodingQueues.begin(),
    qEnd=forwardCodingQueues.end();
  for(; qCur!=qEnd ; ++qCur)
    {
      SignalQueue &queue=**qCur;

      // *********** ITERATE THROUGH THE MAIN QUEUE *************

      SignalQueueIterator sCur=queue.begin(), sEnd=queue.end();
      while(sCur!=sEnd)
	{
	  // Figure out which phase the TAG occurs in, and terminate it
	  Signal &signal=**sCur;
	  int delta=TAGposition-signal.posOfBaseFollowingConsensus();
	  int stoppedPhase=posmod(-delta);
	  Propagator &prop=signal.getPropagator(queue.getContentType());
	  prop[stoppedPhase]=NEGATIVE_INFINITY;

	  // If all phases of this signal are stopped, remove it from queue
	  if(isinf(prop[(stoppedPhase+1)%3]) &&
	     isinf(prop[(stoppedPhase+2)%3]))
	    {
	      SignalQueueIterator victim=sCur;
#ifdef DEBUG
	      cout<<"ECLIPSING "<<signal.getSignalType()<<"@"<<signal.getConsensusPosition()<<endl;
#endif
	      ++sCur; // otherwise queue.erase() will invalidate it
	      queue.drop(victim);
	      //###signal.decReferenceCount(); // for leaving the queue
	      //###if(signal.getReferenceCount()==0) delete &signal;
	      continue;
	    }
	  ++sCur;
	}

      // ********** NOW DO THE SAME THING FOR THE HOLDING QUEUE *********

      sCur=queue.getHoldingQueue();
      sEnd=queue.holdingQueueEnd();
      while(sCur!=sEnd)
	{
	  // Figure out which phase the TAG occurs in, and terminate it
	  Signal &signal=**sCur;
	  int delta=TAGposition-signal.posOfBaseFollowingConsensus();
	  int stoppedPhase=posmod(-delta);
	  Propagator &prop=signal.getPropagator(queue.getContentType());
	  prop[stoppedPhase]=NEGATIVE_INFINITY;

	  // If all phases of this signal are stopped, remove it from queue
	  if(isinf(prop[(stoppedPhase+1)%3]) &&
	     isinf(prop[(stoppedPhase+2)%3]))
	    {
	      SignalQueueIterator victim=sCur;
	      ++sCur; // otherwise queue.erase() will invalidate it
	      queue.dropFromHoldingQueue(victim);
	      //###signal.decReferenceCount(); // for leaving the queue
	      //###if(signal.getReferenceCount()==0) delete &signal;
	      continue;
	    }
	  ++sCur;
	}
    }
}



void GeneZilla::terminateReverseORFs(int TAGposition)
{
  BOOM::Vector<SignalQueue*>::iterator qCur=forwardCodingQueues.begin(),
    qEnd=forwardCodingQueues.end();
  for(; qCur!=qEnd ; ++qCur)
    {
      SignalQueue &queue=**qCur;

      // *********** ITERATE THROUGH THE MAIN QUEUE *************

      SignalQueueIterator sCur=queue.begin(), sEnd=queue.end();
      while(sCur!=sEnd)
	{
	  // Figure out which phase the TAG occurs in, and terminate it
	  Signal &signal=**sCur;
	  int delta=TAGposition-signal.posOfBaseFollowingConsensus();
	  int stoppedPhase=posmod(delta-1);
	  Propagator &prop=signal.getPropagator(queue.getContentType());
	  prop[stoppedPhase]=NEGATIVE_INFINITY;

	  // If all phases of this signal are stopped, remove it from queue
	  if(isinf(prop[(stoppedPhase+1)%3]) &&
	     isinf(prop[(stoppedPhase+2)%3]))
	    {
	      SignalQueueIterator victim=sCur;
	      ++sCur; // otherwise queue.erase() will invalidate it
	      queue.drop(victim);
	      //###signal.decReferenceCount(); // for leaving the queue
	      //###if(signal.getReferenceCount()==0) delete &signal;
	      continue;
	    }
	  ++sCur;
	}

      // ********** NOW DO THE SAME THING FOR THE HOLDING QUEUE *********

      sCur=queue.getHoldingQueue();
      sEnd=queue.holdingQueueEnd();
      while(sCur!=sEnd)
	{
	  // Figure out which phase the TAG occurs in, and terminate it
	  Signal &signal=**sCur;
	  int delta=TAGposition-signal.posOfBaseFollowingConsensus();
	  int stoppedPhase=posmod(delta-1);
	  Propagator &prop=signal.getPropagator(queue.getContentType());
	  prop[stoppedPhase]=NEGATIVE_INFINITY;

	  // If all phases of this signal are stopped, remove it from queue
	  if(isinf(prop[(stoppedPhase+1)%3]) &&
	     isinf(prop[(stoppedPhase+2)%3]))
	    {
	      SignalQueueIterator victim=sCur;
	      ++sCur; // otherwise queue.erase() will invalidate it
	      queue.dropFromHoldingQueue(victim);
	      //###signal.decReferenceCount(); // for leaving the queue
	      //###if(signal.getReferenceCount()==0) delete &signal;
	      continue;
	    }
	  ++sCur;
	}
    }  
}



BOOM::Stack<Signal*> *GeneZilla::traceBack(Signal *rightTerminus,int phase)
{
  cout << "Tracing back to find optimal path..." << endl;

  // Trace back and build stack of visited signals
  BOOM::Stack<Signal*> *stk=new BOOM::Stack<Signal*>;
  stk->push(rightTerminus);
  Signal *signal=rightTerminus;
  while(true)
    {
      //cout<<phase << " " << signal->getSignalType()<<" @ " << signal->getConsensusPosition() << endl;
      Signal *pred=signal->getPredecessor(phase);
      if(!pred) break;
      int distance=signal->getConsensusPosition() -
	pred->posOfBaseFollowingConsensus();
      stk->push(pred);
      switch(pred->getSignalType())
	{
	case NEG_TAG:
	case NEG_POLYA:
	case NEG_ATG:     
	case NEG_PROM:    
	  phase=2;
	  break;
	case PROM:        
	case ATG:
	case TAG:         
	case POLYA:       
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



void GeneZilla::generateGff(BOOM::Stack<Signal*> *path)
{
  cout << "##gff-version 2" << endl;
  cout << "##source-version " << PROGRAM_NAME << " " << VERSION << endl;
  cout << "##date " << getDateAndTime();
  cout << "##Type DNA" << endl;
  //###cout << "# Full-parse score: " << getGeneModelScore(path) << endl << endl;

  // Invariant: the stack always contains at least 2 signals, and both
  // the top signal and the bottom signal are bogus (left & right termini)
  Signal *thisSignal=path->pop();
  int transgrp=0;
  while(!path->isEmpty())
    {
      SignalType thisType=thisSignal->getSignalType();
      Signal *nextSignal=path->pop();
      SignalType nextType=nextSignal->getSignalType();
      int thisPos=thisSignal->getConsensusPosition();
      int nextPos=nextSignal->getConsensusPosition();
#ifdef DEBUG
      cout << thisType << " @ " << thisPos << endl;
#endif

      switch(thisType)
	{
	case ATG: // initial-exon or single-exon
	  ++transgrp;
	  switch(nextType)
	    {
	    case GT:
	      cout << substrateId << '\t' << PROGRAM_NAME << "\tinitial-exon\t"
		   << thisPos+1 << '\t' << nextPos << "\t." << "\t+\t0\t"
		   << "transgrp=" << transgrp << ";" << endl;
	      break;
	    case TAG:
	      cout << substrateId << '\t' << PROGRAM_NAME << "\tsingle-exon\t"
		   << thisPos+1 << '\t' << nextPos+3 << "\t." << "\t+\t0\t"
		   << "transgrp=" << transgrp << ";" << endl;
	      break;
	    }
	  break;
	case AG: // internal-exon or final-exon
	  switch(nextType)
	    {
	    case GT:
	      cout << substrateId << '\t' << PROGRAM_NAME 
		   << "\tinternal-exon\t"
		   << thisPos+3 << '\t' << nextPos << "\t." << "\t+\t0\t"
		   << "transgrp=" << transgrp << ";" << endl;
	      break;
	    case TAG:
	      cout << substrateId << '\t' << PROGRAM_NAME << "\tfinal-exon\t"
		   << thisPos+3 << '\t' << nextPos+3 << "\t." << "\t+\t0\t"
		   << "transgrp=" << transgrp << ";" << endl;
	      break;
	    }
	  break;
	case NEG_TAG: // final-exon or single-exon
	  ++transgrp;
	  switch(nextType)
	    {
	    case NEG_AG:
	      cout << substrateId << '\t' << PROGRAM_NAME << "\tfinal-exon\t"
		   << thisPos+1 << '\t' << nextPos << "\t." << "\t-\t0\t"
		   << "transgrp=" << transgrp << ";" << endl;
	      break;
	    case NEG_ATG:
	      cout << substrateId << '\t' << PROGRAM_NAME << "\tsingle-exon\t"
		   << thisPos+1 << '\t' << nextPos+3 << "\t." << "\t-\t0\t"
		   << "transgrp=" << transgrp << ";" << endl;
	      break;
	    }
	  break;
	case NEG_GT: // internal-exon or initial-exon
	  switch(nextType)
	    {
	    case NEG_AG:
	      cout << substrateId << '\t' << PROGRAM_NAME 
		   << "\tinternal-exon\t"
		   << thisPos+3 << '\t' << nextPos << "\t." << "\t-\t0\t"
		   << "transgrp=" << transgrp << ";" << endl;
	      break;
	    case NEG_ATG:
	      cout << substrateId << '\t' << PROGRAM_NAME << "\tinitial-exon\t"
		   << thisPos+3 << '\t' << nextPos+3 << "\t." << "\t-\t0\t"
		   << "transgrp=" << transgrp << ";" << endl;
	      break;
	    }
	  break;
	}
      thisSignal=nextSignal;
    }
}



