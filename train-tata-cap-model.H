/****************************************************************
 train-tata-cap-model.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_TRAIN_TATA_CAP_MODEL_H
#define INCL_TRAIN_TATA_CAP_MODEL_H

#include <string>
#include <iostream>
#include "BOOM/CommandLine.H"
#include "BOOM/DnaAlphabet.H"
#include "ContentSensor.H"
#include "WMM.H"
#include "GarbageCollector.H"
using namespace std;

class Application
{
  double intergenicModel[5];
  GarbageIgnorer GC;
  WMM *capModel;
  WMM *capIntergenicRatioModel;
  int minSeparation, maxSeparation;

  void loadIntergenic(const BOOM::String &filename);
  void loadCapModel(const BOOM::String &filename);
  void writeOutput(const BOOM::String &tataFile,const BOOM::String &outfile);
public:
  Application();
  int main(int argc,char *argv[]);
};

/*************************************************************************
THE TATA MODEL:
---------------

WMM
PROMOTER
-20.485315	15	5
15	1	4
+
-1.8135644  -0.99730599	-0.92331541	-86.060501	-2.5484858	
-3.0944483  -2.1286318	-3.0386398	-86.060486	-0.23851101	
-0.09684135 -86.060478	-5.1498976	-86.060471	-2.4476109	
-5.1159959  -3.6535122	-5.2400484	-86.060463	-0.037909545	
-0.09310265 -86.060455	-4.2336068	-86.060448	-2.5982993	
-0.36240566 -86.06044	-86.06044	-86.060432	-1.1907276	
-0.08218665 -4.7330036	-2.9374633	-86.060425	-4.0686769	
-0.55077994 -5.5727544	-2.1430206	-86.06041	-1.1960046	
-0.9276045  -2.2036452	-0.89208609	-86.060394	-2.4733734	
-1.9268922  -1.0613165	-0.97312564	-86.060379	-2.036382	
-1.5205116  -0.98993725	-1.1041275	-86.060371	-2.5472076	
-1.5697365  -1.1196316	-1.0972133	-86.060356	-2.0272288	
-1.5389794  -1.2275826	-1.0891571	-86.060341	-1.8585405	
-1.7384083  -1.2884419	-1.021929	-86.060333	-1.6681269	
-1.5925808  -1.3337018	-1.0362011	-86.060318	-1.7242877	


THE CAP-SITE MODEL:
-------------------

WMM
PROMOTER
-10.84373	8	5
8	1	2
+
-1.82201	-1.84263	-1.47972	-inf	-0.79385
-inf	0	-inf	-inf	-inf
-0.05076	-inf	-inf	-inf	-3.00578
-2.45573	-1.31938	-0.96024	-inf	-1.3318
-1.37002	-1.15995	-inf	-inf	-0.83863
-1.50914	-0.94314	-1.88519	-inf	-1.43716
-1.90717	-1.27118	-1.42337	-inf	-1.10866
-1.8018	-1.14948	-1.68848	-inf	-1.09871
**************************************************************************/

#endif
