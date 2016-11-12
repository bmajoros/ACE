/****************************************************************
 LabelMatrix.C
 Copyright (C)2014 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <math.h>
#include <iostream>
#include "LabelMatrix.H"
#include "BOOM/File.H"
using namespace std;
using namespace BOOM;


LabelMatrix::LabelMatrix(const String &filename)
{
  load(filename);
}



float LabelMatrix::operator()(GeneModelLabel from,GeneModelLabel to)
{
  return M[from][to];
}



void LabelMatrix::load(const String &filename)
{
  M.resize(NumGeneModelLabels,NumGeneModelLabels);
  File f(filename);
  f.getline(); // comment line -- ignore
  f.getline(); // header line -- ignore
  for(int i=0 ; i<NumGeneModelLabels ; ++i) {
    String line=f.getline();
    line.trimWhitespace();
    if(line.isEmpty()) continue;
    Vector<String> &fields=*line.getFields();
    if(fields.size()<NumGeneModelLabels+1) 
      throw filename+" : error parsing matrix";
    for(int j=0 ; j<NumGeneModelLabels ; ++j)
      M[i][j]=fields[j+1].asFloat();
    delete &fields;
  }
}



void LabelMatrix::convertToLogs()
{
  const int firstDim=M.getFirstDim(), secondDim=M.getSecondDim();
  for(int x=0 ; x<firstDim ; ++x)
    for(int y=0 ; y<secondDim ; ++y)
      M[x][y]=log(M[x][y]);
}




void LabelMatrix::printOn(ostream &os) const
{
  const int firstDim=M.getFirstDim(), secondDim=M.getSecondDim();
  for(int y=0 ; y<secondDim ; ++y) os<<static_cast<GeneModelLabel>(y)<<"\t";
  os<<endl;
  for(int x=0 ; x<firstDim ; ++x) {
    os<<static_cast<GeneModelLabel>(x)<<"\t";
    for(int y=0 ; y<secondDim ; ++y)
      os<<M[x][y]<<"\t";
    os<<endl;
  }
}



ostream &operator<<(ostream &os,const LabelMatrix &M)
{
  M.printOn(os);
  return os;
}



