/****************************************************************
 GffReader.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "GffReader.H"
//#include <iostream>
//#include "BOOM/File.H"
//#include "genezilla.H"
//#include "GeneZilla.H"


GffReader::GffReader(const BOOM::String &filename,
		     BOOM::Vector<SignalSensor*> &sensors,
		     GarbageCollector &gc)
  : gc(gc)
{
  BOOM::Vector<SignalSensor*>::iterator cur=sensors.begin(), 
    end=sensors.end();
  for(; cur!=end ; ++cur) 
    {
      cout<<"ADDING SENSOR FOR TYPE \""<<(*cur)->getSignalType()<<"\"="<<*cur<<endl;
      signalSensors[(*cur)->getSignalType()]=*cur;
    }

  load(filename);
}



GffReader::~GffReader()
{
  /*
  SignalList::iterator lCur=signals.begin(), lEnd=signals.end();
  for(; lCur!=lEnd ; ++lCur) delete *lCur;
  */

  BOOM::Map<int,SignalList*>::iterator hCur=posToSignals.begin(),
    hEnd=posToSignals.end();
  for(; hCur!=hEnd ; ++hCur) delete (*hCur).second;
}



GffReader::SignalList &GffReader::getSignals()
{
  return signals;
}



GffReader::SignalList *GffReader::signalsAt(int position)
{
  if(!posToSignals.isDefined(position)) return NULL;
  return posToSignals[position];
}



void GffReader::load(const BOOM::String &filename)
{
  BOOM::File file(filename);
  while(!file.eof())
    {
      BOOM::String line=file.readLine();
      if(file.eof()) break;
      BOOM::Vector<BOOM::String> &fields=*line.getFields();
      int numFields=fields.size();
      if(numFields>=9 && fields[0][0]!='#')
	{
	  BOOM::String &featureType=fields[2];
	  int begin=fields[3].asInt()-1;
	  int end=fields[4].asInt();
	  float score=fields[5].asFloat();
	  Strand strand=fields[6][0];
	  int contextWindowPos;
	  SignalSensor *sensor;
	  switch(strand)
	    {
	    case FORWARD_STRAND:
	      if(featureType=="initial-exon")
		{
		  // ATG
		  sensor=signalSensors[ATG];
		  cout<<"1sensor="<<sensor<<endl;
		  contextWindowPos=begin-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);

		  // GT
		  sensor=signalSensors[GT];
		  cout<<"2sensor="<<sensor<<endl;
		  contextWindowPos=end-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		}
	      else if(featureType=="internal-exon")
		{
		  // AG
		  sensor=signalSensors[AG];
		  cout<<"3sensor="<<sensor<<endl;
		  contextWindowPos=begin-2-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		  
		  // GT
		  sensor=signalSensors[GT];
		  cout<<"4sensor="<<sensor<<endl;
		  contextWindowPos=end-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		}
	      else if(featureType=="final-exon")
		{
		  // AG
		  sensor=signalSensors[AG];
		  cout<<"5sensor="<<sensor<<endl;
		  contextWindowPos=begin-2-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		  
		  // TAG
		  sensor=signalSensors[TAG];
		  cout<<"6sensor="<<sensor<<endl;
		  contextWindowPos=end-3-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		}
	      else if(featureType=="single-exon")
		{
		  // ATG
		  sensor=signalSensors[ATG];
		  cout<<"7sensor="<<sensor<<endl;
		  contextWindowPos=begin-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		  
		  // TAG
		  sensor=signalSensors[TAG];
		  cout<<"8sensor="<<sensor<<endl;
		  contextWindowPos=end-3-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		}
	      break;

	    //===============================================================
	    case REVERSE_STRAND:

	      if(featureType=="initial-exon")
		{
		  // NEG_GT
		  sensor=signalSensors[NEG_GT];
		  cout<<"9sensor="<<sensor<<endl;
		  contextWindowPos=begin-2-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);

		  // NEG_ATG
		  sensor=signalSensors[NEG_ATG];
		  cout<<"10sensor="<<sensor<<endl;
		  contextWindowPos=end-3-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		}
	      else if(featureType=="internal-exon")
		{
		  // NEG_GT
		  sensor=signalSensors[NEG_GT];
		  cout<<"11sensor="<<sensor<<endl;
		  contextWindowPos=begin-2-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);

		  // NEG_AG
		  sensor=signalSensors[NEG_AG];
		  cout<<"12sensor="<<sensor<<endl;
		  contextWindowPos=end-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		}
	      else if(featureType=="final-exon")
		{
		  // NEG_TAG
		  sensor=signalSensors[NEG_TAG];
		  cout<<"13sensor="<<sensor<<endl;
		  contextWindowPos=begin-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);

		  // NEG_AG
		  sensor=signalSensors[NEG_AG];
		  cout<<"14sensor="<<sensor<<endl;
		  contextWindowPos=end-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		}
	      else if(featureType=="single-exon")
		{
		  // NEG_TAG
		  sensor=signalSensors[NEG_TAG];
		  cout<<"15sensor="<<sensor<<endl;
		  contextWindowPos=begin-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);

		  // NEG_ATG
		  sensor=signalSensors[NEG_ATG];
		  cout<<"16sensor="<<sensor<<endl;
		  contextWindowPos=end-3-sensor->getConsensusOffset();
		  newSignal(contextWindowPos,score,*sensor);
		}
	      break;
	    }
	}
      delete &fields;
    }
}



void GffReader::newSignal(int contextWindowPos,float score,
			  SignalSensor &sensor)
{
  SignalPtr signal=new Signal(contextWindowPos,score,sensor,gc);
  linkBack(signal);
  signals.push_back(signal);
  if(!posToSignals.isDefined(contextWindowPos))
    posToSignals[contextWindowPos]=new SignalList;
  posToSignals[contextWindowPos]->push_back(signal);
}


/*
void GffReader::createLeftTerminus(SignalType signalType)
{
  switch(signalType)
    {
    case ATG:
    case TAG:
    case GT:
    case AG:
    case PROM:
    case POLYA:
    case NEG_ATG:
    case NEG_TAG:
    case NEG_GT:
    case NEG_AG:
    case NEG_PROM:
    case NEG_POLYA:
      
    }
}
*/



int GffReader::getPhase(SignalPtr pred)
{
  // ### THIS IS ONLY VALID IF NO PARTIAL GENES ARE PROCESSED...

  int phase=(pred->getStrand()==FORWARD_STRAND ? 0 : 2);
  for(int i=0 ; i<3 ; ++i)
    if(pred->getPredecessor(i))
      phase=i;
  return phase;
}



void GffReader::linkBack(SignalPtr signal)
{
  if(signals.isEmpty()) return;
  int pos=signal->getConsensusPosition();
  SignalPtr pred=signals.getLast();
  int predPos=pred->posOfBaseFollowingConsensus();
  int predPhase=getPhase(pred);
  int phase;
  switch(signal->getSignalType())
    {
      // + strand:
    case ATG:         phase=0;break;
    case TAG:         phase=0;break;
    case PROM:        phase=0;break;
    case POLYA:       phase=0;break;
    case GT:          phase=(predPhase+pos-predPos)%3;break;
    case AG:          phase=predPhase;break;

      // - strand:
    case NEG_ATG:     phase=2;break;
    case NEG_TAG:     phase=2;break;
    case NEG_PROM:    phase=2;break;
    case NEG_POLYA:   phase=2;break;
    case NEG_GT:      phase=predPhase;break;
    case NEG_AG:      phase=posmod(predPhase-pos+predPos);break;
    }

#ifdef DEBUG
  //cout<<"GENSCAN----> PHASE OF "<<signal.getSignalType()<<" @ "<<signal.getConsensusPosition()<<" = "<<phase<<endl;
#endif

  signal->setPredecessor(phase,pred);
}



SignalPtr GffReader::getSignal(int position,SignalType t)
{
  SignalPtr pred;
  return getSignal(position,t,pred);
}



SignalPtr GffReader::getSignal(int position,SignalType t,SignalPtr &pred)
{
  pred=NULL;
  SignalSensor *sensor=signalSensors[t];
  int wndPos=position-sensor->getConsensusOffset();
  SignalList *sigs=signalsAt(wndPos);
  if(!sigs) return NULL;
  int n=sigs->size();
  for(int i=0 ; i<n ; ++i)
    {
      SignalPtr genscanSig=(*sigs)[i];
      if(genscanSig->getSignalType()==t)
	{
	  int gsPhase=getPhase(genscanSig);
	  pred=genscanSig->getPredecessor(gsPhase);
	  return genscanSig;
	}
    }
  return NULL;
}



void GffReader::addLeftTerminus(SignalPtr s)
{
  signals.push_front(s);
  int contextWindowPos=s->getContextWindowPosition();
  if(!posToSignals.isDefined(contextWindowPos))
    posToSignals[contextWindowPos]=new SignalList;
  posToSignals[contextWindowPos]->push_back(s);
}



void GffReader::addRightTerminus(SignalPtr s)
{
  signals.push_back(s);
  int contextWindowPos=s->getContextWindowPosition();
  if(!posToSignals.isDefined(contextWindowPos))
    posToSignals[contextWindowPos]=new SignalList;
  posToSignals[contextWindowPos]->push_back(s);
}



#ifdef EXPLICIT_GRAPHS
ParseGraph *GffReader::toGraph(GeneZilla &tigrScan)
{
  ParseGraph *graph=new ParseGraph(tigrScan);
  BOOM::Vector<SignalPtr>::iterator cur=signals.begin(), end=signals.end();
  for(; cur!=end ; ++cur)
    {
      SignalPtr signal=*cur;
      graph->addVertex(signal);
    }
  signals.clear();
  return graph;
}
#endif

