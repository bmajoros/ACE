/****************************************************************
 Transitions.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "Transitions.H"
#include <math.h>
#include <iostream>
#include <fstream>
#include "BOOM/Constants.H"
#include "genezilla.H"


Transitions::Transitions(int numSignalTypes,istream &is,float optimism,
			 float intronOptimism)
  : matrix(numSignalTypes,numSignalTypes), optimism(optimism),
    intronOptimism(intronOptimism)
{
  matrix.setAllTo(0);
  load(is);
}



void Transitions::load(istream &is)
{
  while(!is.eof()) {
    BOOM::String line;
    getline(is,line);
    line.trimWhitespace();
    if(!line.empty() && line[0]=='#') continue;
    BOOM::Vector<BOOM::String> fields;
    line.getFields(fields," \t");
    if(fields.size()==5) {
      SignalType from=stringToSignalType(fields[0]);
      SignalType to=stringToSignalType(fields[2]);
      double P=fields[4].asDouble();
      matrix[from][to]=P;
    }
    else if(fields.size()>0)
      throw BOOM::String("syntax error in transition file: ")+line;
  }
  //cout<<matrix<<endl;
  normalize();
  //cout<<"normalized:\n"<<matrix<<endl;
  logify();
  //cout<<"logged:\n"<<matrix<<endl;
}



void Transitions::normalize()
{
  const int N=matrix.getFirstDim();
  for(int x=0 ; x<N ; ++x) {
    float sum=0;
    for(int y=0 ; y<N ; ++y) sum+=matrix[x][y];
    for(int y=0 ; y<N ; ++y) matrix[x][y]/=sum;
  }
}



void Transitions::logify()
{
  const int N=matrix.getFirstDim();
  for(int x=0 ; x<N ; ++x)
    for(int y=0 ; y<N ; ++y) 
      matrix[x][y]=log(matrix[x][y]);
}



double Transitions::getLogP(SignalType from,SignalType to)
{
  if(beginsCoding(to)) 
    return 
      optimism +
      matrix[static_cast<int>(from)][static_cast<int>(to)];
  else if(beginsIntron(to))
    return 
      intronOptimism +
      matrix[static_cast<int>(from)][static_cast<int>(to)];    
  else
    return 
      matrix[static_cast<int>(from)][static_cast<int>(to)];
}

