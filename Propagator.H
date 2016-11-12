/****************************************************************
 Propagator.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/

#ifndef INCL_Propagator_H
#define INCL_Propagator_H

#include <iostream>
#include "BOOM/Strand.H"
#include "Signal.H"
#include "BOOM/Constants.H"
class SignalQueue;

/***************************************************************
 PROPAGATOR SEMANTICS

 A propagator consists of an array of three values and a current
 position.  It is used to track the inductive score at individual
 base positions of the path passing through a given signal.  If 
 signal S has propagator P which is currently at position N, then 
 P[X] denotes the inductive score at position N given that the
 base at position N lies on a path which passes through signal S
 in phase X.  Thus, P[] can be used to determine which phase is
 best for signal S.  Just remember that for P[X], X is the phase
 of the signal S, not of the base at position N!
 
***************************************************************/

class Propagator
{
  friend class SignalQueue;
  int position; // necessary because signal lengths vary
protected:
  double scores[3]; // initialized to contentScore (all 3 phases)
public:
  inline Propagator(int position);
  virtual ~Propagator() {}
  inline double &operator[](int);
  inline int getPosition();
  inline void incrementPosition();
  inline void update(Propagator &accumulator,Strand,bool isCoding,
		     SignalPtr mySignal);
  inline void resetScoresForCoding(); // sets all 3 scores to 0
  inline void resetScoresForNoncoding(); // [0]=0, [1]=-inf, [2]=-inf
  inline void addFrameZeroAndAdvance(double d) {*scores+=d;++position;}
  void printOn(ostream &);
  void resetPosition() {position=-1;}
};

ostream &operator<<(ostream &,Propagator &);


inline Propagator::Propagator(int position)
  : position(position)
{
}




inline double &Propagator::operator[](int i)
{
  return scores[i];
}



inline int Propagator::getPosition()
{
  return position;
}



inline void Propagator::incrementPosition()
{
  ++position;
}



inline void Propagator::resetScoresForCoding() 
{
  double *sc=scores;
  *sc=*(sc+1)=*(sc+2)=0;    // scores[0 & 1 & 2]
}



inline void Propagator::resetScoresForNoncoding()
{
  double *sc=scores;
  *sc=0;                                 // scores[0]
  *(sc+1)=*(sc+2)=NEGATIVE_INFINITY;     // scores[1 & 2]
}



inline void Propagator::update(Propagator &accum,Strand strand,
			       bool isCoding,SignalPtr mySignal)
{
  //cout<<"Propagator::update() isCoding="<<isCoding<<endl;
  if(isCoding)
    {
#ifdef DEBUG
	  if(isinf(accum.scores[0]) || isinf(accum.scores[1]) || isinf(accum.scores[2])) cout<<"XXXXXXXXX ADDING -INF TO CODING PROPAGATOR"<<endl;
#endif

      int frame=mySignal->frameOfBaseFollowingConsensus();
      double *sc=scores;
      switch(strand)
	{
	case FORWARD_STRAND:
	  // ### THIS CAN BE OPTIMIZED FURTHER BY COMPUTING POSMOD(0-FRAME)
	  //   AND THEN DUPLICATING THE CODE FOR THE THREE POSSIBLE RESULTS,
	  //   THEREBY ALLOWING A POINTER TO BE USED FOR ACCUM.SCORES

	  //cout<<"updating scores, frame="<<frame<<" 0-frame="<<posmod(0-frame)<<" 1-frame="<<posmod(1-frame)<<" 2-frame="<<posmod(2-frame)<<endl;
	  //cout<<*sc<<" "<<*(sc+1)<<" "<<*(sc+2)<<" => ";
	  *sc+=accum.scores[posmod(0-frame)];     // scores[0]
	  *(sc+1)+=accum.scores[posmod(1-frame)]; // scores[1]
	  *(sc+2)+=accum.scores[posmod(2-frame)]; // scores[2]
	  //cout<<*sc<<" "<<*(sc+1)<<" "<<*(sc+2)<<endl;
	  //cout<<"accum was: "<<accum.scores[0]<<" "<<accum.scores[1]<<" "<<accum.scores[2]<<endl;
	  /*
	  {
	    int i=posmod(-frame);
	    switch(i)
	      {
	      case 0:
		{
		  double *as=accum.scores;
		  *sc+=*(as++);     // scores[0]
		  *(sc+1)+=*(as++); // scores[1]
		  *(sc+2)+=*as; // scores[2]
		}
		break;
	      case 1:
		{
		  double *as=accum.scores;
		  *(sc+2)+=*(as++); // scores[2]
		  *sc+=*(as++);     // scores[0]
		  *(sc+1)+=*as; // scores[1]
		}
		break;
	      case 2:
		{
		  double *as=accum.scores;
		  *(sc+1)+=*(as++); // scores[1]
		  *(sc+2)+=*(as++); // scores[2]
		  *sc+=*as;     // scores[0]
		}
		break;
	      }
	  }
	  */
	  break;
	case REVERSE_STRAND:
	  // ### optimize this just like on the forward strand
	  *sc+=accum.scores[posmod(0+frame)];     // scores[0]
	  *(sc+1)+=accum.scores[posmod(1+frame)]; // scores[1]
	  *(sc+2)+=accum.scores[posmod(2+frame)]; // scores[2]
	  break;
	}
    }
  else // NONCODING
    {
      // Copy accumulator slot 0 into all three slots of this propagator
      double accumScore0=accum.scores[0];
      double *sc=scores;
      *sc+=accumScore0;
      *(sc+1)+=accumScore0;
      *(sc+2)+=accumScore0;  
    }
  position=accum.position;
}



inline ostream &operator<<(ostream &os,Propagator &prop)
{
  prop.printOn(os);
  return os;
}



inline void Propagator::printOn(ostream &os)
{
  os << position << ":" << "[" << scores[0]
     << "," << scores[1] << "," << scores[2]
     << "]";
}




#endif
