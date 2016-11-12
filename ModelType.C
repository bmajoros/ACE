/*
 ModelType.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include "ModelType.H"
#include "BOOM/Map.H"
#include "BOOM/String.H"

static BOOM::Map<BOOM::String,ModelType> modelTypeMap;
ModelTypeInitializer ModelTypeInitializer::initializer;


ModelTypeInitializer::ModelTypeInitializer()
{
  modelTypeMap["WMM"]=WMM_MODEL;
  modelTypeMap["WAM"]=WAM_MODEL;
  modelTypeMap["WWAM"]=WWAM_MODEL;
  modelTypeMap["3P"]=THREE_PERIODIC;
  modelTypeMap["MC"]=MARKOV_CHAIN;
  modelTypeMap["FastMC"]=MARKOV_CHAIN;
  modelTypeMap["Fast3PMC"]=MARKOV_CHAIN;
  modelTypeMap["IMM"]=IMM_MODEL;
  modelTypeMap["CB"]=CODON_BIAS;
  modelTypeMap["MDD"]=MDD;
  modelTypeMap["HMM"]=HMM_MODEL;
  modelTypeMap["NSMC"]=NONSTATIONARY_MC;
  modelTypeMap["3PIMM"]=IMM_3P;
  modelTypeMap["BranchAcceptor"]=BRANCH_ACCEPTOR;
}



BOOM::String toString(ModelType t)
{
  switch(t)
    {
    case WMM_MODEL: return "WMM";
    case WAM_MODEL: return "WAM";
    case WWAM_MODEL: return "WWAM"; 
    case THREE_PERIODIC: return "3P";
    case FAST_MC: return "FastMC";
    case FAST_3PMC: return "Fast3PMC";
    case MARKOV_CHAIN: return "MC";
    case IMM_MODEL: return "IMM";
    case CODON_BIAS: return "CB";
    case MDD: return "MDD";
    case SIGNAL_MODEL: return "SIGNAL_MODEL";
    case HMM_MODEL: return "HMM";
    case NONSTATIONARY_MC: return "NSMC";
    case IMM_3P: return "3PIMM";
    case BRANCH_ACCEPTOR: return "BranchAcceptor";
    }
}



ModelType stringToModelType(const BOOM::String &s)
{
  if(!modelTypeMap.isDefined(s)) 
    throw BOOM::String("undefined ModelType: ")+s;
  return modelTypeMap[s];
}

