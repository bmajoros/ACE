#include "BOOM/SumLogProbs.H"
using namespace std;
using namespace BOOM;

void test()
{
  Vector<double> v;
  v.push_back(1e-5);
  v.push_back(1.001e-5);
  v.push_back(1.03e-5);
  v.push_back(1.01e-5);

  // Compute sums in log and non-log space
  double sum=0.0;
  Vector<double> logs;
  for(Vector<double>::iterator cur=v.begin(), end=v.end() ; cur!=end ; ++cur) {
    sum+=*cur;
    logs.push_back(log(*cur));
  }
  double sumLog=sumLogProbs(logs);
  cout<<"sum="<<sum<<" sum from logs="<<exp(sumLog)<<endl;
  
  //###
  logs.clear();
  logs.push_back(-10800185);
  logs.push_back(-10800186);
  //logs.push_back(-10800187);
  //logs.push_back(-10800188);
  sumLog=sumLogProbs(logs);
  sumLog=sumLogProbs(logs[0],logs[1]);
  sumLog=logs[0]+log(1+exp(logs[1]-logs[0]));
  cout.precision(12);
  cout<<logs[1]-logs[0]<<endl;
  cout<<exp(logs[1]-logs[0])<<endl;
  cout<<1+exp(logs[1]-logs[0])<<endl;
  cout<<log(1+exp(logs[1]-logs[0]))<<endl;
  cout<<"sumLog="<<sumLog<<endl;
  //###

  // Convert to probabilities
  /*
  float rawSum=0.0;
  for(Vector<float>::iterator cur=v.begin(), end=v.end() ; cur!=end ; ++cur) {
    float p=*cur;
    p/=sum;
    rawSum+=p;
  }  
  cout<<"raw sum="<<rawSum<<endl;
  */
  double sumFromLogs=0.0;
  for(Vector<double>::iterator cur=logs.begin(), end=logs.end() ; 
      cur!=end ; ++cur) {
    double logP=*cur;
    double p=exp(logP-sumLog);
    cout<<"p="<<p<<endl;
    sumFromLogs+=p;
  }
  cout<<"sum from logs="<<sumFromLogs<<endl;
}


int main(int argc,char *argv[])
{
  try {
    test();
    return 0;
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    { cerr << "STL exception caught in main:\n" << e.what() << endl; }
  catch(...)
    { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}


