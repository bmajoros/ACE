/**************************************************************
 Partition.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#ifndef INCL_Partition_H
#define INCL_Partition_H
#include <iostream>
#include "BOOM/Set.H"
#include "BOOM/Vector.H"
#include "TrainingSequence.H"

enum Direction
  {
    DIR_LEFT,
    DIR_RIGHT
  };

class Partition
{
  int index;
  BOOM::Set<Symbol> leftResidues;

  void load(istream &);
public:
  Partition();
  Partition(Partition &,int index);
  Partition(istream &);
  double split(BOOM::Vector<TrainingSequence*> &parent,
	       BOOM::Vector<TrainingSequence*> &leftChild,
	       BOOM::Vector<TrainingSequence*> &rightChild);
  void addSymbol(Symbol);
  void setIndex(int);
  int getIndex();
  bool isLeftSymbol(Symbol);
  void printOn(ostream &) const;
  Direction getDirection(const Sequence &,int windowBegin);
  void save(ostream &);
  Partition *reverseComplement(int sequenceLength);
  Partition *clone() const;
};

ostream &operator<<(ostream &,const Partition &);

#endif
