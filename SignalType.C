/*
 SignalType.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 */
#include "SignalType.H"
#include "BOOM/Map.H"
#include "BOOM/String.H"
#include "BOOM/Exceptions.H"

static BOOM::Map<BOOM::String,SignalType> typeNames;
static BOOM::Map<BOOM::String,SignalType> longTypeNames;
TypeNamesInitializer TypeNamesInitializer::initializer;


TypeNamesInitializer::TypeNamesInitializer()
{
  typeNames["LT"]=LEFT_TERMINUS;
  typeNames["RT"]=RIGHT_TERMINUS;
  typeNames["ATG"]=ATG;
  typeNames["TAG"]=TAG;
  typeNames["GT"]=GT;
  typeNames["AG"]=AG;
  typeNames["UTR5GT"]=UTR5GT;
  typeNames["UTR5AG"]=UTR5AG;
  typeNames["UTR3GT"]=UTR3GT;
  typeNames["UTR3AG"]=UTR3AG;
  typeNames["-UTR5GT"]=NEG_UTR5GT;
  typeNames["-UTR5AG"]=NEG_UTR5AG;
  typeNames["-UTR3GT"]=NEG_UTR3GT;
  typeNames["-UTR3AG"]=NEG_UTR3AG;
  typeNames["TSS"]=TSS;
  typeNames["TES"]=TES;
  typeNames["-ATG"]=NEG_ATG;
  typeNames["-TAG"]=NEG_TAG;
  typeNames["-GT"]=NEG_GT;
  typeNames["-AG"]=NEG_AG;
  typeNames["-TSS"]=NEG_TSS;
  typeNames["-TES"]=NEG_TES;
  typeNames["GT0"]=GT0;
  typeNames["GT1"]=GT1;
  typeNames["GT2"]=GT2;
  typeNames["AG0"]=AG0;
  typeNames["AG1"]=AG1;
  typeNames["AG2"]=AG2;
  typeNames["-GT0"]=NEG_GT0;
  typeNames["-GT1"]=NEG_GT1;
  typeNames["-GT2"]=NEG_GT2;
  typeNames["-AG0"]=NEG_AG0;
  typeNames["-AG1"]=NEG_AG1;
  typeNames["-AG2"]=NEG_AG2;

  longTypeNames["start-codon"]=ATG;
  longTypeNames["stop-codon"]=TAG;
  longTypeNames["donor"]=GT;
  longTypeNames["acceptor"]=AG;
  longTypeNames["transcription-start-site"]=TSS;
  longTypeNames["transcription-end-site"]=TES;
  longTypeNames["left-terminus"]=LEFT_TERMINUS;
  longTypeNames["right-terminus"]=RIGHT_TERMINUS;
}



SignalType stringToSignalType(const BOOM::String &s)
{
  if(typeNames.isDefined(s)) return typeNames[s];
  if(longTypeNames.isDefined(s)) return longTypeNames[s];
  throw BOOM::String("undefined SignalType: ")+s;
}



BOOM::String signalTypeToString(SignalType t)
{
  switch(t)
    {
    case LEFT_TERMINUS:  return "LT";
    case RIGHT_TERMINUS: return "RT";
    case ATG:            return "ATG";
    case TAG:            return "TAG";
    case GT:             return "GT"; 
    case GT0:             return "GT0"; 
    case GT1:             return "GT1"; 
    case GT2:             return "GT2"; 
    case AG:             return "AG";
    case AG0:             return "AG0";
    case AG1:             return "AG1";
    case AG2:             return "AG2";
    case UTR5GT:         return "UTR5GT";
    case UTR5AG:         return "UTR5AG";
    case UTR3GT:         return "UTR3GT";
    case UTR3AG:         return "UTR3AG";
    case TSS:            return "TSS";
    case TES:            return "TES";
    case NEG_ATG:        return "-ATG";    
    case NEG_TAG:        return "-TAG";    
    case NEG_GT:         return "-GT";     
    case NEG_GT0:         return "-GT0";     
    case NEG_GT1:         return "-GT1";     
    case NEG_GT2:         return "-GT2";     
    case NEG_AG:         return "-AG";     
    case NEG_AG0:         return "-AG0";     
    case NEG_AG1:         return "-AG1";     
    case NEG_AG2:         return "-AG2";     
    case NEG_TSS:        return "-TSS";
    case NEG_TES:        return "-TES";  
    case NEG_UTR5GT:     return "-UTR5GT";
    case NEG_UTR5AG:     return "-UTR5AG";
    case NEG_UTR3GT:     return "-UTR3GT";
    case NEG_UTR3AG:     return "-UTR3AG";
    }
  throw "signalTypeToString()";
}



BOOM::String signalTypeToName(SignalType t)
{
  switch(t)
    {
    case LEFT_TERMINUS:  return "left-terminus";
    case RIGHT_TERMINUS: return "right-terminus";
    case ATG:            return "start-codon";
    case TAG:            return "stop-codon";
    case GT0: case GT1: case GT2: case UTR5GT: case UTR3GT:
    case GT:             return "donor"; 
    case AG0: case AG1: case AG2: case UTR5AG: case UTR3AG:
    case AG:             return "acceptor"; 
      //case UTR5GT:         return "UTR5-donor";
      //case UTR5AG:         return "UTR5-acceptor";
      //case UTR3GT:         return "UTR3-donor";
      //case UTR3AG:         return "UTR3-acceptor";
    case TSS:            return "transcription-start-site";
    case TES:            return "transcription-end-site";   
    case NEG_ATG:        return "start-codon";    
    case NEG_TAG:        return "stop-codon";    
    case NEG_GT0: case NEG_GT1: case NEG_GT2:
    case NEG_GT:         return "donor";     
    case NEG_AG0: case NEG_AG1: case NEG_AG2:
    case NEG_AG:         return "acceptor";     
    case NEG_TSS:        return "transcription-start-site";
    case NEG_TES:        return "transcription-end-site";
    case NEG_UTR5GT:     return "UTR5-donor";
    case NEG_UTR5AG:     return "UTR5-acceptor";
    case NEG_UTR3GT:     return "UTR3-donor";
    case NEG_UTR3AG:     return "UTR3-acceptor";
    }
  throw "signalTypeToName()";
}



ostream &operator<<(ostream &os,SignalType t)
{
  switch(t)
    {
    case LEFT_TERMINUS:os<< "LT";           break;
    case RIGHT_TERMINUS:os<<"RT";           break;
    case ATG:         os << "ATG";          break;
    case TAG:         os << "TAG";          break;
    case GT:          os << "GT";           break;
    case GT0:         os << "GT0";          break;
    case GT1:         os << "GT1";          break;
    case GT2:         os << "GT2";          break;
    case AG:          os << "AG";           break;
    case AG0:         os << "AG0";          break;
    case AG1:         os << "AG1";          break;
    case AG2:         os << "AG2";          break;
    case UTR5GT:      os << "UTR5GT";       break;
    case UTR5AG:      os << "UTR5AG";       break;
    case UTR3GT:      os << "UTR3GT";       break;
    case UTR3AG:      os << "UTR3AG";       break;
    case TSS:         os << "TSS";          break;
    case TES:         os << "TES";          break;
    case NEG_ATG:     os << "-ATG";         break;
    case NEG_TAG:     os << "-TAG";         break;
    case NEG_GT:      os << "-GT";          break;
    case NEG_GT0:     os << "-GT0";         break;
    case NEG_GT1:     os << "-GT1";         break;
    case NEG_GT2:     os << "-GT2";         break;
    case NEG_AG:      os << "-AG";          break;
    case NEG_AG0:     os << "-AG0";         break;
    case NEG_AG1:     os << "-AG1";         break;
    case NEG_AG2:     os << "-AG2";         break;
    case NEG_TSS:     os << "-TSS";         break;
    case NEG_TES:     os << "-TES";         break;
    case NEG_UTR5GT:  os << "-UTR5GT";      break;
    case NEG_UTR5AG:  os << "-UTR5AG";      break;
    case NEG_UTR3GT:  os << "-UTR3GT";      break;
    case NEG_UTR3AG:  os << "-UTR3AG";      break;
    default: throw BOOM::String("INVALID SIGNAL TYPE: ")+int(t);
    }
  return os;
}



istream &operator>>(istream &is,SignalType &t)
{
  BOOM::String buf;
  is >> buf;
  if(!typeNames.isDefined(buf)) 
    throw BOOM::String(" undefined SignalType: ")+buf;
  t=typeNames[buf];
  return is;
}



SignalType reverseComplement(SignalType t)
{
  switch(t)
    {
    case LEFT_TERMINUS:   return LEFT_TERMINUS;  
    case RIGHT_TERMINUS:  return RIGHT_TERMINUS;
    case ATG:             return NEG_ATG;
    case TAG:             return NEG_TAG;
    case GT:              return NEG_GT;
    case GT0:             return NEG_GT0;
    case GT1:             return NEG_GT1;
    case GT2:             return NEG_GT2;
    case AG:              return NEG_AG;
    case AG0:             return NEG_AG0;
    case AG1:             return NEG_AG1;
    case AG2:             return NEG_AG2;
    case UTR5GT:          return NEG_UTR5GT;
    case UTR5AG:          return NEG_UTR5AG;
    case UTR3GT:          return NEG_UTR3GT;
    case UTR3AG:          return NEG_UTR3AG;
    case TSS:             return NEG_TSS;
    case TES:             return NEG_TES;
    case NEG_ATG:         return ATG;
    case NEG_TAG:         return TAG;
    case NEG_GT:          return GT;
    case NEG_AG:          return AG;
    case NEG_TSS:         return TSS;
    case NEG_TES:         return TES;
    case NEG_UTR5GT:      return UTR5GT;
    case NEG_UTR5AG:      return UTR5AG;
    case NEG_UTR3GT:      return UTR3GT;
    case NEG_UTR3AG:      return UTR3AG;
    }
}



Strand getStrand(SignalType t)
{
  switch(t)
    {
    case NO_SIGNAL_TYPE:
    case LEFT_TERMINUS:
    case RIGHT_TERMINUS:
      return NO_STRAND;

    case ATG:
    case TAG:
    case GT:
    case AG:
    case UTR5GT:
    case UTR5AG:
    case UTR3GT:
    case UTR3AG:
    case TSS:
    case TES:
    case GT0:
    case GT1:
    case GT2:
    case AG0:
    case AG1:
    case AG2:
        return FORWARD_STRAND;

    case NEG_ATG:
    case NEG_TAG:
    case NEG_GT:
    case NEG_AG:
    case NEG_UTR5GT:
    case NEG_UTR5AG:
    case NEG_UTR3GT:
    case NEG_UTR3AG:
    case NEG_TSS:
    case NEG_TES:
    case NEG_GT0:
    case NEG_GT1:
    case NEG_GT2:
    case NEG_AG0:
    case NEG_AG1:
    case NEG_AG2:
      return REVERSE_STRAND;
    }
}



bool endsCoding(SignalType t)
{
  switch(t)
    {
    case TAG:
    case GT: case GT0: case GT1: case GT2:
    case UTR5GT:
    case UTR3GT:
    case NEG_ATG:
    case NEG_AG: case NEG_AG0: case NEG_AG1: case NEG_AG2:
    case NEG_UTR5AG:
    case NEG_UTR3AG:
      return true;
    }
  return false;
}


bool beginsCoding(SignalType t)
{
  switch(t)
    {
    case ATG: 
    case AG: case AG0: case AG1: case AG2:
    case UTR5AG:
    case UTR3AG:
    case NEG_TAG:
    case NEG_GT: case NEG_GT0: case NEG_GT1: case NEG_GT2:
    case NEG_UTR5GT:
    case NEG_UTR3GT:
      return true;
    }
  return false;
}


bool beginsIntron(SignalType t)
{
  switch(t)
    {
    case GT: case GT0: case GT1: case GT2:
    case UTR5GT:
    case UTR3GT:
    case NEG_AG: case NEG_AG0: case NEG_AG1: case NEG_AG2:
    case NEG_UTR5AG:
    case NEG_UTR3AG:
      return true;
    }
  return false;
}


bool endsIntron(SignalType t)
{
  switch(t)
    {
    case NEG_GT: case NEG_GT0: case NEG_GT1: case NEG_GT2:
    case NEG_UTR5GT:
    case NEG_UTR3GT:
    case AG: case AG0: case AG1: case AG2:
    case UTR5AG:
    case UTR3AG:
      return true;
    }
  return false;
}



SignalType dropStrand(SignalType t)
{
  switch(t)
    {
    case LEFT_TERMINUS:
    case RIGHT_TERMINUS:
    case ATG:
    case TAG:
    case GT: case GT0: case GT1: case GT2:
    case AG: case AG0: case AG1: case AG2:
    case UTR5GT:
    case UTR5AG:
    case UTR3GT:
    case UTR3AG:
    case TSS:
    case TES:
      return t;
    case NEG_ATG:
    case NEG_TAG:
    case NEG_GT: case NEG_GT0: case NEG_GT1: case NEG_GT2:
    case NEG_AG: case NEG_AG0: case NEG_AG1: case NEG_AG2:
    case NEG_UTR5GT:
    case NEG_UTR5AG:
    case NEG_UTR3GT:
    case NEG_UTR3AG:
    case NEG_TSS:
    case NEG_TES:
      return reverseComplement(t);
    }
  INTERNAL_ERROR;
}



