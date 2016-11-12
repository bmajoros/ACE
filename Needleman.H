/****************************************************************
 Needleman.H
 bmajoros@duke.edu

 Implements an affine version of the the Needleman-Wunsch global 
 alignment algorithm.  Can perform "partial" alignment -- i.e., 
 where terminal indels are not penalized.
 ****************************************************************/
#ifndef INCL_Needleman_H
#define INCL_Needleman_H
#include "AlignmentPath.H"
#include "BOOM/Array2D.H"
#include "AlignmentSubstMatrix.H"

template<class FLT_TYP>
class Needleman {
public:
  Needleman(Alphabet &,const Sequence &,const Sequence &,
		const AlignmentSubstMatrix<FLT_TYP> &,
		FLT_TYP gapOpenPenalty,FLT_TYP gapExtendPenalty);
  AlignmentPath *fullAlignment(bool ungapped=false);
  AlignmentPath *partialAlignment();
  FLT_TYP getPartialScore();
  enum Direction {LEFT,DOWN,DIAGONAL};
private:
  struct Cell {
    FLT_TYP scores[3];  
    Direction links[3];
    Direction bestDir();
    Cell() { scores[0]=scores[1]=scores[2]=NEGATIVE_INFINITY; }
  };
  Alphabet &alphabet;
  const Sequence &s1, &s2;
  Array2D<Cell> alignmentMatrix; // index like this: M[column][row]
  const AlignmentSubstMatrix<FLT_TYP> &substitutionMatrix;
  const FLT_TYP gapOpenPenalty, gapExtendPenalty;
  AlignmentPath *getAlignment(int finalRow,int finalColumn);
  AlignmentPath *align(FLT_TYP finalOpenPenalty,
		   FLT_TYP finalExtendPenalty,bool ungapped=false);
  AlignmentPath *align_ungapped(FLT_TYP finalOpenPenalty,
                   FLT_TYP finalExtendPenalty);                                           
  void linkBack(Cell &thisCell,Cell &leftCell,Cell &downCell,
		Cell &diagonalCell,FLT_TYP gapOpenPenalty,
		FLT_TYP gapExtendPenalty,int x,int y);
  void linkBack(Cell &thisCell,                                                             
                Cell &diagonalCell,FLT_TYP gapOpenPenalty,
                FLT_TYP gapExtendPenalty,int x,int y);
  void followLink(Direction,int &x,int &y);
};

template<class FLT_TYP>
ostream &operator<<(ostream &,typename Needleman<FLT_TYP>::Direction);



#include "BOOM/Stack.H"
#include "BOOM/Constants.H"


template<class FLT_TYP>
Needleman<FLT_TYP>::Needleman(Alphabet &alphabet,
				      const Sequence &s1,
				      const Sequence &s2,
				      const AlignmentSubstMatrix<FLT_TYP> &subst,
				      FLT_TYP gapOpenPenalty,
				      FLT_TYP gapExtendPenalty)
  : s1(s1), s2(s2), alphabet(alphabet), 
    alignmentMatrix(s1.getLength()+1,s2.getLength()+1),
    gapOpenPenalty(gapOpenPenalty), gapExtendPenalty(gapExtendPenalty),
    substitutionMatrix(subst)
{
  if(s1.getLength()<1 || s2.getLength()<1)
    throw "Aligning to a zero-length sequence in Needleman::Needleman";
}



template<class FLT_TYP>
typename Needleman<FLT_TYP>::Direction 
  Needleman<FLT_TYP>::Cell::bestDir()
{
  FLT_TYP left=scores[Needleman<FLT_TYP>::LEFT]; 
  FLT_TYP down=scores[Needleman<FLT_TYP>::DOWN];
  FLT_TYP diag=scores[Needleman<FLT_TYP>::DIAGONAL];
  if(left>down)
    if(left>diag)
      return Needleman<FLT_TYP>::LEFT;
    else
      return Needleman<FLT_TYP>::DIAGONAL;
  else
    if(down>diag)
      return Needleman<FLT_TYP>::DOWN;
    else
      return Needleman<FLT_TYP>::DIAGONAL;
}



template<class FLT_TYP>
ostream &operator<<(ostream &os,
		    typename Needleman<FLT_TYP>::Direction dir)
{
  switch(dir)
    {
    case Needleman<FLT_TYP>::LEFT:
      os<<"LEFT"<<endl;
      break;
    case Needleman<FLT_TYP>::DOWN:
      os<<"DOWN"<<endl;
      break;
    case Needleman<FLT_TYP>::DIAGONAL:
      os<<"DIAGONAL"<<endl;
      break;
    }
}



template<class FLT_TYP>
AlignmentPath *Needleman<FLT_TYP>::getAlignment(int x,int y)
{
  // Trace back from the end cell along the optimal path to
  // the origin in the lower left corner of the DP matrix
  Cell &cell=alignmentMatrix[x][y];
  Direction dir=cell.bestDir();
  FLT_TYP alignmentScore=cell.scores[dir];

  Stack<Direction> pathStack;
  while(x!=0 || y!=0) {
    Cell &cell=alignmentMatrix[x][y];
    pathStack.push(dir);
    followLink(dir,x,y);
    dir=cell.links[dir];
  }

  // Reverse the stack contents to produce a forward-oriented
  // alignment
  AlignmentPath &alignment=*new AlignmentPath(s1,s2,alphabet,alignmentScore);
  const int len=pathStack.size();
  for(int i=0 ; i<len ; ++i) {
    Direction dir=pathStack.pop();
    switch(dir)
      {
      case LEFT:
	alignment+=FIRST_UNMATCHED;
	break;
      case DOWN:
	alignment+=SECOND_UNMATCHED;
	break;
      case DIAGONAL:
	alignment+=MATCH;
	break;
      }
  }
  return &alignment;
}



template<class FLT_TYP>
AlignmentPath *Needleman<FLT_TYP>::fullAlignment(bool ungapped)
{
  return align(gapOpenPenalty,gapExtendPenalty,ungapped);
}



template<class FLT_TYP>
AlignmentPath *Needleman<FLT_TYP>::partialAlignment()
{
  return align(0,0);
}



template<class FLT_TYP>
AlignmentPath *Needleman<FLT_TYP>::align_ungapped(FLT_TYP finalOpenPenalty,
						  FLT_TYP finalExtendPenalty)
{
  int matrixWidth=alignmentMatrix.getFirstDim();
  int matrixHeight=alignmentMatrix.getSecondDim();

  // Initialize the bottom row and leftmost column of the dynamic
  // programming matrix
  Cell &origin=alignmentMatrix[0][0];
  for(int i=0 ; i<3 ; ++i) origin.scores[i]=0;
  for(int i=1 ; i<matrixWidth ; ++i) {
    Cell &bottomCell=alignmentMatrix[i][0];
    bottomCell.links[LEFT]=LEFT;
    bottomCell.scores[LEFT]=
      alignmentMatrix[i-1][0].scores[LEFT]+gapExtendPenalty;
    bottomCell.scores[DOWN]=bottomCell.scores[DIAGONAL]=NEGATIVE_INFINITY;
  }
  for(int i=1 ; i<matrixHeight ; ++i) {
    Cell &leftCell=alignmentMatrix[0][i];
    leftCell.links[DOWN]=DOWN;
    leftCell.scores[DOWN]=
      alignmentMatrix[0][i-1].scores[DOWN]+gapExtendPenalty;
    leftCell.scores[LEFT]=leftCell.scores[DIAGONAL]=NEGATIVE_INFINITY;
  }

  // Fill out the main portion of the matrix
  int widthMinusOne=matrixWidth-1, heightMinusOne=matrixHeight-1;
  for(int x=1 ; x<widthMinusOne ; ++x)
    for(int y=1 ; y<heightMinusOne ; ++y) {
      Cell &thisCell=alignmentMatrix[x][y],
	&diagonalCell=alignmentMatrix[x-1][y-1];
      linkBack(thisCell,diagonalCell,
	       gapOpenPenalty,gapExtendPenalty,x,y);
    }

  // Finish up by filling out the top row and rightmost column
  // (these use a *potentially* different set of gap penalties to
  // accommodate "partial" alignments)
  for(int i=1 ; i<matrixWidth ; ++i) {
    Cell &topCell=alignmentMatrix[i][heightMinusOne],
      topCellLeft=alignmentMatrix[i-1][heightMinusOne],
      topCellDown=alignmentMatrix[i][heightMinusOne-1],
      topCellDiagonal=alignmentMatrix[i-1][heightMinusOne-1];
    linkBack(topCell,
	     topCellLeft,topCellDown,topCellDiagonal,
	     finalOpenPenalty,finalExtendPenalty,i,heightMinusOne);
  }
  for(int i=1 ; i<matrixHeight ; ++i) {
    Cell &rightmostCell=alignmentMatrix[widthMinusOne][i],
      &rightmostCellLeft=alignmentMatrix[widthMinusOne-1][i],
      &rightmostCellDown=alignmentMatrix[widthMinusOne][i-1],
      &rightmostCellDiagonal=alignmentMatrix[widthMinusOne-1][i-1];
    linkBack(rightmostCell,
	     rightmostCellLeft,rightmostCellDown,
	     rightmostCellDiagonal,finalOpenPenalty,finalExtendPenalty,
	     widthMinusOne,i);
  }
  return getAlignment(widthMinusOne,heightMinusOne);
}
 


template<class FLT_TYP>
AlignmentPath *Needleman<FLT_TYP>::align(FLT_TYP finalOpenPenalty,
					 FLT_TYP finalExtendPenalty,
					 bool ungapped)
{
  if(ungapped) return align_ungapped(finalOpenPenalty,finalExtendPenalty);                                   
  int matrixWidth=alignmentMatrix.getFirstDim();
  int matrixHeight=alignmentMatrix.getSecondDim();

  // Initialize the bottom row and leftmost column of the dynamic 
  // programming matrix
  Cell &origin=alignmentMatrix[0][0];
  for(int i=0 ; i<3 ; ++i) origin.scores[i]=0;
  for(int i=1 ; i<matrixWidth ; ++i) {
    Cell &bottomCell=alignmentMatrix[i][0];
    bottomCell.links[LEFT]=LEFT;
    bottomCell.scores[LEFT]=
      alignmentMatrix[i-1][0].scores[LEFT]+finalExtendPenalty;//gapExtendPenalty;
    bottomCell.scores[DOWN]=bottomCell.scores[DIAGONAL]=NEGATIVE_INFINITY;
  }
  for(int i=1 ; i<matrixHeight ; ++i) {
    Cell &leftCell=alignmentMatrix[0][i];
    leftCell.links[DOWN]=DOWN;
    leftCell.scores[DOWN]=
      alignmentMatrix[0][i-1].scores[DOWN]+finalExtendPenalty;//gapExtendPenalty;
    leftCell.scores[LEFT]=leftCell.scores[DIAGONAL]=NEGATIVE_INFINITY;
  }

  // Fill out the main portion of the matrix
  int widthMinusOne=matrixWidth-1, heightMinusOne=matrixHeight-1;
  for(int x=1 ; x<widthMinusOne ; ++x)
    for(int y=1 ; y<heightMinusOne ; ++y) {
      Cell &thisCell=alignmentMatrix[x][y],
	&leftCell=alignmentMatrix[x-1][y],
	&downCell=alignmentMatrix[x][y-1],
	&diagonalCell=alignmentMatrix[x-1][y-1];
      linkBack(thisCell,leftCell,downCell,diagonalCell,
	       gapOpenPenalty,gapExtendPenalty,x,y);
    }

  // Finish up by filling out the top row and rightmost column
  // (these use a *potentially* different set of gap penalties to
  // accommodate "partial" alignments)
  for(int i=1 ; i<matrixWidth ; ++i) {
    Cell &topCell=alignmentMatrix[i][heightMinusOne],
      topCellLeft=alignmentMatrix[i-1][heightMinusOne],
      topCellDown=alignmentMatrix[i][heightMinusOne-1],
      topCellDiagonal=alignmentMatrix[i-1][heightMinusOne-1];
    linkBack(topCell,topCellLeft,topCellDown,topCellDiagonal,
	     finalOpenPenalty,finalExtendPenalty,i,heightMinusOne);
  }
  for(int i=1 ; i<matrixHeight ; ++i) {
    Cell &rightmostCell=alignmentMatrix[widthMinusOne][i],
      &rightmostCellLeft=alignmentMatrix[widthMinusOne-1][i],
      &rightmostCellDown=alignmentMatrix[widthMinusOne][i-1],
      &rightmostCellDiagonal=alignmentMatrix[widthMinusOne-1][i-1];
    linkBack(rightmostCell,rightmostCellLeft,rightmostCellDown,
	     rightmostCellDiagonal,finalOpenPenalty,finalExtendPenalty,
	     widthMinusOne,i);
  }

  return getAlignment(widthMinusOne,heightMinusOne);
}



/************************************************************************
 This method handles the actual linking back of one cell to a predecessor
 cell, either immediately left, below, or diagonally of this cell.  Note
 that we do not allow insertions directly next to deletions -- at least
 one match must come between adjacent insertions and deletions, although
 this could could easily be changed to relax this restriction.
 */
template<class FLT_TYP>
void Needleman<FLT_TYP>::linkBack(Cell &thisCell,Cell &leftCell,
				  Cell &downCell,
				  Cell &diagonalCell,
				  FLT_TYP gapOpenPenalty,
				  FLT_TYP gapExtendPenalty,
				  int x,int y)
{
  // INSERTION
  FLT_TYP leftLeft=leftCell.scores[LEFT]+gapExtendPenalty;
  FLT_TYP leftDiag=leftCell.scores[DIAGONAL]+gapOpenPenalty;
  if(leftLeft>leftDiag) {
    thisCell.scores[LEFT]=leftLeft;
    thisCell.links[LEFT]=LEFT;
  }
  else {
    thisCell.scores[LEFT]=leftDiag;
    thisCell.links[LEFT]=DIAGONAL;
  }

  // DOWN
  FLT_TYP downDown=downCell.scores[DOWN]+gapExtendPenalty;
  FLT_TYP downDiag=downCell.scores[DIAGONAL]+gapOpenPenalty;
  if(downDown>downDiag) {
    thisCell.scores[DOWN]=downDown;
    thisCell.links[DOWN]=DOWN;
  }
  else {
    thisCell.scores[DOWN]=downDiag;
    thisCell.links[DOWN]=DIAGONAL;
  }
  
  // DIAGONAL
  FLT_TYP matchScore=substitutionMatrix(s1[x-1],s2[y-1]);
  FLT_TYP diagLeft=diagonalCell.scores[LEFT]+matchScore;
  FLT_TYP diagDown=diagonalCell.scores[DOWN]+matchScore;
  FLT_TYP diagDiag=diagonalCell.scores[DIAGONAL]+matchScore;
  if(diagLeft>diagDown)
    if(diagLeft>diagDiag) { // left is best
      thisCell.scores[DIAGONAL]=diagLeft;
      thisCell.links[DIAGONAL]=LEFT;
    }
    else goto DIAGDIAG;
  else 
    if(diagDown>diagDiag) { // down is best
      thisCell.scores[DIAGONAL]=diagDown;
      thisCell.links[DIAGONAL]=DOWN;
    }
    else { // diagonal is best
    DIAGDIAG:
      thisCell.scores[DIAGONAL]=diagDiag;
      thisCell.links[DIAGONAL]=DIAGONAL;
    }
}



template<class FLT_TYP>
void Needleman<FLT_TYP>::linkBack(Cell &thisCell,
				  Cell &diagonalCell,
				  FLT_TYP gapOpenPenalty,
				  FLT_TYP gapExtendPenalty,
				  int x,int y)
{
  thisCell.scores[LEFT]=thisCell.scores[DOWN]=NEGATIVE_INFINITY;
  thisCell.links[LEFT]=thisCell.links[DOWN]=DIAGONAL;

  // DIAGONAL
  FLT_TYP matchScore=substitutionMatrix(s1[x-1],s2[y-1]);
  FLT_TYP diagLeft=diagonalCell.scores[LEFT]+matchScore;
  FLT_TYP diagDown=diagonalCell.scores[DOWN]+matchScore;
  FLT_TYP diagDiag=diagonalCell.scores[DIAGONAL]+matchScore;
  if(diagLeft>diagDown)
    if(diagLeft>diagDiag) { // left is best
      thisCell.scores[DIAGONAL]=diagLeft;
      thisCell.links[DIAGONAL]=LEFT;
    }
    else goto DIAGDIAG;
  else
    if(diagDown>diagDiag) { // down is best
      thisCell.scores[DIAGONAL]=diagDown;
      thisCell.links[DIAGONAL]=DOWN;
    }
    else { // diagonal is best
     DIAGDIAG:
      thisCell.scores[DIAGONAL]=diagDiag;
      thisCell.links[DIAGONAL]=DIAGONAL;
    }
} 



template<class FLT_TYP>
void Needleman<FLT_TYP>::followLink(Direction dir,int &x,int &y)
{
  switch(dir)
    {
    case LEFT:
      --x;
      break;
    case DOWN:
      --y;
      break;
    case DIAGONAL:
      --x;
      --y;
      break;
    }
}



/**********************************************************************
 Get the score of an alignment ending at the last diagonal along the
 optimal path.  To do this we start at the upper, right cell and trace
 back until we get to the first diagonal.  In this way, we don't penalize
 for indels at the end of the alignment.
 */
template<class FLT_TYP>
FLT_TYP Needleman<FLT_TYP>::getPartialScore()
{
  int x=alignmentMatrix.getFirstDim()-1;
  int y=alignmentMatrix.getSecondDim()-1;
  while(x!=0 || y!=0) {
    Cell &cell=alignmentMatrix[x][y];
    Direction dir=cell.bestDir();
    if(dir==DIAGONAL) return cell.scores[dir];
    followLink(dir,x,y);
  }
  return 0;
}


#endif
