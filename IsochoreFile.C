/****************************************************************
 IsochoreFile.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "IsochoreFile.H"
#include <iostream>
#include "BOOM/File.H"

// G+C <= 100% : genezilla.cfg
BOOM::Regex IsochoreFile::lineRegex(
  "G[+]C\\s*<=\\s*(\\S+)%\\s*:\\s*(\\S+)");


BOOM::ConfigFile *IsochoreFile::parse(const BOOM::String &filename,
				    double gcContent,
				    BOOM::String &gffMessage)
{
  gcContent*=100;
  BOOM::File infile(filename);
  while(!infile.eof())
    {
      BOOM::String line=infile.readLine();
      if(lineRegex.search(line))
	{
	  double upperBound=lineRegex[1].asDouble();
	  BOOM::String configFilename=lineRegex[2];
	  if(gcContent<=upperBound)
	    {
	      cerr << "\tG+C%=" << gcContent << " <= " << upperBound
		   << " -> using " << configFilename << endl;
	      gffMessage=BOOM::String("# G+C=")+int(10*gcContent)/10.0+
		"% -- used isochore file: "+configFilename;
	      return new BOOM::ConfigFile(configFilename);
	    }
	}
    }
  throw BOOM::String("Error parsing ")+filename+" : G+C content = "+
    gcContent+"%";
}




