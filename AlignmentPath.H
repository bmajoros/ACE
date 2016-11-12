/****************************************************************
 AlignmentPath.H
 bmajoros@duke.edu

 Represents an alignment between two sequences (nucleotide, amino
 acid, or other).  The alignment is represented by a list of
 match classifications (MATCH/INSERTION/DELETION) that describe
 the individual points on the path through an alignment matrix 
 (i.e., DIAGONAL/UP/RIGHT).  This is a (relatively) compact 
 representation that is bounded by the sum of the two sequence 
 lengths, but it does not give constant time access to a random 
 base.  It is time-efficient only when processing the entire 
 alignment from left to right.  A constant-time random-access 
 alignment representation could be derived from this.
 ****************************************************************/
#ifndef INCL_AlignmentPath_H
#define INCL_AlignmentPath_H
#include "BOOM/Sequence.H"
#include <iostream>
#include "MatchType.H"
#include "AlignmentSubstMatrix.H"
using namespace std;
using namespace BOOM;

class AlignmentPath {
public:
  static int MAX_WIDTH;
  AlignmentPath(const Sequence &,const Sequence &,Alphabet &,double score=0);
  const Sequence &getFirstSequence() const;
  const Sequence &getSecondSequence() const;
  int getAlignmentLength() const; // length of *full* alignment
  int getAlignedLength() const; // not including indels at either end
  MatchType operator[](int position) const;
  AlignmentPath &operator+=(MatchType);
  double getScore() const;
  void countMismatches(int &mismatches,int &insertions) const;
  int countNearMatches(AlignmentSubstMatrix<float> &);
  int countNearMatches(AlignmentSubstMatrix<double> &);
  void getResidualsOnRight(Sequence &,Sequence &);
  void countNonNColumns(int &cols,int &matches) const;
  void getMatchExtent(int &firstBegin,int &firstEnd,int &secondBegin,int &secondEnd) const;
  String getCigarString() const;
private:
  const Sequence &s1, &s2;
  Vector<MatchType> matchData;
  double score;
  Alphabet &alphabet;
  friend ostream &operator<<(ostream &,const AlignmentPath &);
  void printOn(ostream &) const;
};

ostream &operator<<(ostream &,const AlignmentPath &);

#endif
