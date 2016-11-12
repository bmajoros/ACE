/**************************************************************
 Partition.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
***************************************************************/
#include "Partition.H"
#include <iostream>
#include "BOOM/Alphabet.H"

extern Alphabet alphabet;

Partition::Partition(istream &is)
  : index(-1)
{
  // ctor

  load(is);
}



Partition::Partition()
  : index(-1)
{
  // ctor
}



Partition::Partition(Partition &other,int index)
  : index(index), leftResidues(other.leftResidues)
{
  // ctor
}



double Partition::split(BOOM::Vector<TrainingSequence*> &parent,
			BOOM::Vector<TrainingSequence*> &leftChild,
			BOOM::Vector<TrainingSequence*> &rightChild)
{
  int n=parent.size();
  for(int i=0 ; i<n ; ++i)
    {
      TrainingSequence &seq=*parent[i];
      if(leftResidues.isMember(seq[index]))
	leftChild.push_back(&seq);
      else
	rightChild.push_back(&seq);
    }
  return leftChild.size()/double(n);
}



Direction Partition::getDirection(const Sequence &seq,int begin)
{
  int trueIndex=begin+index;
  //###if(trueIndex>=seq.getLength()) return DIR_LEFT;

  return leftResidues.isMember(seq[trueIndex]) ? DIR_LEFT : DIR_RIGHT;
}



void Partition::addSymbol(Symbol s)
{
  leftResidues.insert(s);
}



void Partition::setIndex(int i)
{
  index=i;
}



int Partition::getIndex()
{
  return index;
}



bool Partition::isLeftSymbol(Symbol s)
{
  return leftResidues.isMember(s);
}



void Partition::printOn(ostream &os) const
{
  os << "#" << index << " in {";
  BOOM::Set<Symbol>::iterator cur=leftResidues.begin(),
    end=leftResidues.end();
  for(; cur!=end ; ++cur)
    {
      Symbol s=*cur;
      os << alphabet.lookup(s);
    }
  os << "}";
}



ostream &operator<<(ostream &os,const Partition &partition)
{
  partition.printOn(os);
  return os;
}


void Partition::save(ostream &os)
{
  os << index << "\t" << leftResidues.size() << "\t";
  BOOM::Set<Symbol>::iterator cur=leftResidues.begin(),
    end=leftResidues.end();
  for(; cur!=end ; ++cur)
    {
      Symbol s=*cur;
      os << alphabet.lookup(s);
    }
  os << endl;
}



void Partition::load(istream &is)
{
  int numResidues;
  is >> index >> numResidues;
  for(int i=0 ; i<numResidues ; ++i)
    {
      char c;
      is >> c;
      if(isspace(c)) {--i;continue;}
      Symbol s=alphabet.lookup(c);
      leftResidues.insert(s);
    }
}



Partition *Partition::reverseComplement(int sequenceLength)
{
  Partition *newPartition=new Partition;
  newPartition->index=sequenceLength-1-index;
  BOOM::Set<Symbol> &newSet=newPartition->leftResidues;
  BOOM::Set<Symbol>::iterator cur=leftResidues.begin(), end=leftResidues.end();
  for(; cur!=end ; ++cur)
    newSet.insert(alphabet.complement(*cur));
  return newPartition;
}



Partition *Partition::clone() const
{
  return new Partition(*this);
}

