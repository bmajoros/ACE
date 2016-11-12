/****************************************************************
 train-tata-cap-model.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include "train-tata-cap-model.H"
#include "BOOM/Constants.H"

Alphabet alphabet;


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    BOOM::CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=6) 
      throw BOOM::String(
"\n\n\
train-tata-cap-model <tata.model> <cap.model> <intergenic.model> \n\
                      <min-submodel-separation> <max-submodel-separation>\n\
                      <outfile>\n\n");
    BOOM::String tataFile=cmd.arg(0);
    BOOM::String capFile=cmd.arg(1);
    BOOM::String intergenicFile=cmd.arg(2);
    minSeparation=cmd.arg(3).asInt();
    maxSeparation=cmd.arg(4).asInt();
    BOOM::String outfile=cmd.arg(5);

    // Load the intergenic model and reduce to zeroth order
    loadIntergenic(intergenicFile);

    // Load the cap site model and produce a likelihood ratio version of it
    loadCapModel(capFile);

    // Write out the submodels into a single model file
    writeOutput(tataFile,outfile);
    
    return 0;
  }




void Application::loadIntergenic(const BOOM::String &filename)
{
  Alphabet &alphabet=DnaAlphabet::global();
  ContentSensor *fullModel=ContentSensor::load(filename);
  char base[2];
  base[1]='\0';

  for(int s=0 ; s<5 ; ++s)
    {
      base[0]=alphabet.lookup(s);
      Sequence seq(base,alphabet);
      intergenicModel[s]=
	fullModel->scoreSingleBase(seq,BOOM::String(base),0,seq[0],base[0]);
    }
  delete fullModel;
}



void Application::loadCapModel(const BOOM::String &filename)
{
  // Load the cap model
  capModel=dynamic_cast<WMM*>(SignalSensor::load(filename,GC));

  // Compute a version of the cap model which uses likelihood ratios
  // (cap / intergenic)
  capIntergenicRatioModel=new WMM(GC,*capModel);
  BOOM::Array2D<float> &matrix=capIntergenicRatioModel->getMatrix();
  int windowSize=matrix.getFirstDim();
  for(int pos=0 ; pos<windowSize ; ++pos)
    {
      BOOM::Array2D<float>::RowIn2DArray<float> column=matrix[pos];
      for(int s=0 ; s<5 ; ++s)
	column[s]-=intergenicModel[s];
    }
}



void Application::writeOutput(const BOOM::String &tataFile,
			      const BOOM::String &outfile)
{
  // Load TATA model
  SignalSensor *tata=SignalSensor::load(tataFile,GC);

  // Create output file and write header
  ofstream os(outfile.c_str());
  os.precision(8);
  os<<"TataCapModel"<<endl;
  os<<minSeparation<<"\t"<<maxSeparation<<endl;
  
  // Write out the TATA model
  tata->save(os);

  // Write out the intergenic model
  Alphabet &alphabet=DnaAlphabet::global();
  os<<"MC\nINTERGENIC\n0\t0\t1\n5"<<endl;
  os<<"A\n"<<intergenicModel[alphabet.lookup('A')]<<endl;
  os<<"C\n"<<intergenicModel[alphabet.lookup('C')]<<endl;
  os<<"G\n"<<intergenicModel[alphabet.lookup('G')]<<endl;
  os<<"N\n"<<intergenicModel[alphabet.lookup('N')]<<endl;
  os<<"T\n"<<intergenicModel[alphabet.lookup('T')]<<endl;
  os<<"MC\nINTERGENIC\n0\t0\t1\n5"<<endl;
  os<<"A\n"<<intergenicModel[alphabet.lookup('T')]<<endl;
  os<<"C\n"<<intergenicModel[alphabet.lookup('G')]<<endl;
  os<<"G\n"<<intergenicModel[alphabet.lookup('C')]<<endl;
  os<<"N\n"<<intergenicModel[alphabet.lookup('N')]<<endl;
  os<<"T\n"<<intergenicModel[alphabet.lookup('A')]<<endl;

  // Write out the CAP model
  capModel->save(os);

  // Write out the CAP/intergenic ratio model
  capIntergenicRatioModel->save(os);
}



