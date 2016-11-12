/****************************************************************
 AlignmentPath.C
 bmajoros@duke.edu
 ****************************************************************/
#include "AlignmentPath.H"
#include "BOOM/Exceptions.H"
#include <iostream>
using namespace std;
using namespace BOOM;

int AlignmentPath::MAX_WIDTH=60;



AlignmentPath::AlignmentPath(const Sequence &s1,const Sequence &s2,Alphabet &alpha,
		     double score)
  : s1(s1), s2(s2), alphabet(alpha), score(score)
{
}



const Sequence &AlignmentPath::getFirstSequence() const
{
  return s1;
}



const Sequence &AlignmentPath::getSecondSequence() const
{
  return s2;
}



int AlignmentPath::getAlignmentLength() const
{
  return matchData.size();
}


int AlignmentPath::getAlignedLength() const
{
  
  // Find first match position
  int i, L=matchData.size();
  for(i=0 ; i<L ; ++i) if(matchData[i]==MATCH) break;
  if(i==L) return 0; // no matches!
  int firstMatch=i;

  // Find last match position
  for(i=L-1 ; i>0 ; --i) if(matchData[i]==MATCH) break;
  int lastMatch=i;

  // Compute length of aligned region between first and last match
  return lastMatch-firstMatch+1;
}



MatchType AlignmentPath::operator[](int position) const
{
  return matchData[position];
}



AlignmentPath &AlignmentPath::operator+=(MatchType m)
{
  matchData.push_back(m);
  return *this;
}



void AlignmentPath::printOn(ostream &os) const
{
  String topRow, middleRow, bottomRow;

  Vector<MatchType>::const_iterator cur=matchData.begin(), 
    end=matchData.end();
  int index1=0, index2=0;
  for(; cur!=end ; ++cur)
    switch(*cur)
      {
      case MATCH:
	{
	  Symbol symbol1=s1[index1], symbol2=s2[index2];
	  topRow+=alphabet.lookup(symbol1);
	  bottomRow+=alphabet.lookup(symbol2);
	  middleRow+=(symbol1==symbol2 ? '|' : '*');
	  ++index1;
	  ++index2;
	}
	break;
      case FIRST_UNMATCHED:
	topRow+=alphabet.lookup(s1[index1]);
	bottomRow+='-';
	middleRow+=' ';
	++index1;
	break;
      case SECOND_UNMATCHED:
	topRow+='-';
	bottomRow+=alphabet.lookup(s2[index2]);
	middleRow+=' ';
	++index2;
	break;
      }

  int len=getAlignmentLength(), begin=0, subLen=MAX_WIDTH;
  while(begin<len)
    {
      if(begin+subLen>len) subLen=len-begin;
      os << "Query: " << topRow.substring(begin,subLen) << endl;
      os << "       " << middleRow.substring(begin,subLen) << endl;
      os << "Sbjct: " << bottomRow.substring(begin,subLen) << endl;
      begin+=subLen;
      if(begin<len) os << endl;
    }
}



ostream &operator<<(ostream &os,const AlignmentPath &alignment)
{
  alignment.printOn(os);
  return os;
}



double AlignmentPath::getScore() const
{
  return score;
}


void AlignmentPath::countNonNColumns(int &cols,int &matches) const {
  cols=matches=0;
  Symbol N=alphabet.lookup('N');
  Vector<MatchType>::const_iterator cur=matchData.begin(),
    end=matchData.end();
  int index1=0, index2=0;
  for(; cur!=end ; ++cur) {
    switch(*cur)
      {
      case MATCH:
        {
          Symbol symbol1=s1[index1], symbol2=s2[index2];
	  if(!(symbol1==N) && !(symbol2==N)) {
	    ++cols;
	    if(symbol1==symbol2) ++matches;
	  }
          ++index1;
          ++index2;
        }
        break;
      case FIRST_UNMATCHED:
        ++index1;
        break;
      case SECOND_UNMATCHED:
        ++index2;
        break;
      }
  }
}

void AlignmentPath::countMismatches(int &mismatches,int &insertions) const
{
  mismatches=insertions=0;
  Vector<MatchType>::const_iterator cur=matchData.begin(), 
    end=matchData.end();
  int index1=0, index2=0;
  for(; cur!=end ; ++cur)
    switch(*cur)
      {
      case MATCH:
	{
	  Symbol symbol1=s1[index1], symbol2=s2[index2];
	  if(!(symbol1==symbol2)) ++mismatches;
	  ++index1;
	  ++index2;
	}
	break;
      case FIRST_UNMATCHED:
	++insertions;
	++index1;
	break;
      case SECOND_UNMATCHED:
	++index2;
	++insertions;
	break;
      }
}



int AlignmentPath::countNearMatches(AlignmentSubstMatrix<float> &M)
{
  int count=0;
  Vector<MatchType>::const_iterator cur=matchData.begin(), 
    end=matchData.end();
  int index1=0, index2=0;
  for(; cur!=end ; ++cur)
    switch(*cur)
      {
      case MATCH:
	{
	  Symbol symbol1=s1[index1], symbol2=s2[index2];
	  if(M(symbol1,symbol2)>0) ++count;
	  ++index1;
	  ++index2;
	}
	break;
      case FIRST_UNMATCHED:
	++index1;
	break;
      case SECOND_UNMATCHED:
	++index2;
	break;
      }
  return count;
}



int AlignmentPath::countNearMatches(AlignmentSubstMatrix<double> &M)
{
  int count=0;
  Vector<MatchType>::const_iterator cur=matchData.begin(), 
    end=matchData.end();
  int index1=0, index2=0;
  for(; cur!=end ; ++cur)
    switch(*cur)
      {
      case MATCH:
	{
	  Symbol symbol1=s1[index1], symbol2=s2[index2];
	  if(M(symbol1,symbol2)>0) ++count;
	  ++index1;
	  ++index2;
	}
	break;
      case FIRST_UNMATCHED:
	++index1;
	break;
      case SECOND_UNMATCHED:
	++index2;
	break;
      }
  return count;
}



void AlignmentPath::getResidualsOnRight(Sequence &seq1,Sequence &seq2)
{
  seq1.clear();
  seq2.clear();
  int len=matchData.size();
  int index1=s1.getLength()-1, index2=s2.getLength()-1;
  for(int i=len-1 ; i>0 ; --i)
    switch(matchData[i])
      {
      case MATCH:
	return;
      case FIRST_UNMATCHED:
	seq1.prepend(s1[index1]);
	--index1;
	break;
      case SECOND_UNMATCHED:
	seq2.prepend(s2[index2]);
	--index2;
      }
}


void AlignmentPath::getMatchExtent(int &firstBegin,int &firstEnd,int &secondBegin,int &secondEnd) const
{
  int L=matchData.size();
  int firstPos=0, secondPos=0;
  firstBegin=secondBegin=firstEnd=secondEnd=-1;
  for(int i=0 ; i<L ; ++i) {
    MatchType m=matchData[i];
    switch(m) {
    case MATCH: 
      if(firstBegin<0) firstBegin=firstPos;
      if(secondBegin<0) secondBegin=secondPos;
      firstEnd=firstPos;
      secondEnd=secondPos;
      ++firstPos; ++secondPos; 
      break;
    case FIRST_UNMATCHED: 
      ++firstPos; 
      break;
    case SECOND_UNMATCHED: 
      ++secondPos; 
      break;
    default: INTERNAL_ERROR;
    }
  }
}



char matchTypeToCigarLetter(MatchType m)
{
  switch(m)
    {
    case MATCH: return 'M';
    case FIRST_UNMATCHED: return 'D';
    case SECOND_UNMATCHED: return 'I';
    default: INTERNAL_ERROR;
    }
}



String AlignmentPath::getCigarString() const
{
  int L=matchData.size();
  if(L==0) return "";
  String cigar;
  int opLen=1;
  MatchType matchType=matchData[0];
  char op=matchTypeToCigarLetter(matchType);
  for(int i=1 ; i<L ; ) {
    while(matchData[i]==matchType) { ++opLen; ++i; }
    cigar+=opLen; cigar+=op;
    if(i<L) { 
      matchType=matchData[i]; 
      op=matchTypeToCigarLetter(matchType); 
      opLen=1; }
  }
  return cigar;
}

