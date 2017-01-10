/****************************************************************
 ACEplus.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "ACEplus.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 main()
 ****************************************************************/
int main(int argc,char *argv[])
{
  try {
    ACEplus ace;
    return ace.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    { cerr << "STL exception caught in main:\n" << e.what() << endl; }
  catch(...)
    { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}


