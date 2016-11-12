/****************************************************************
 genezilla.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/


/******************************************************************
  Leave this commented out unless you have a set of signals that
  you definitely want GeneZilla to find (even if the signal sensor
  gives them a low score):
*/
//#define FORCE_SPECIFIC_SIGNALS


/******************************************************************
  Comment or uncomment the following line as per the needs of
  individual projects.  Once you change the status of this line,
  however, you *MUST* issue a "make clean" and recompile the entire 
  project!!! 
*/
//#define EXPLICIT_GRAPHS


/****************************************************************
  Comment out the following line if you do not want to see memory
  usage statistics when the program runs:
*/
#define REPORT_MEMORY_USAGE


/****************************************************************
  Comment out the following line if you do not want to see regular
  progress messages as the program runs:
*/
#define REPORT_PROGRESS



/****************************************************************
  Uncomment the following line if you want to see debugging
  messages while the program is running:
*/
//#define DEBUG


/****************************************************************
  Comment the following line if you want assert() statements to
  be activated (i.e., if you are debugging):
*/
#define NDEBUG



#ifndef INCL_genezilla_H
#define INCL_genezilla_H

// If employing graphs, then use garbage collection, otherwise smart pointers:
#ifdef EXPLICIT_GRAPHS
#undef USE_SMART_POINTERS
#else
#define USE_SMART_POINTERS
//#error BAD!
#endif

#include "BOOM/String.H"
#include "BOOM/Sequence.H"
#include "BOOM/Set.H"
#include "BOOM/SmartPointer.H"
#include "SignalType.H"
#include "ContentType.H"

extern Alphabet alphabet;

class Signal;
#ifndef SignalPtr
typedef BOOM::SmartPointer<Signal> SignalPtr;
#endif

inline int posmod(int); // positive modulus (because % in C++ can be negative)
int mapPhaseBack(int phase,SignalPtr right,SignalPtr left);

int posmod(int x)
{
  int f=x%3;
  if(f>=0) return f;
  return f+3;
}


#endif
