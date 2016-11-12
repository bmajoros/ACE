/****************************************************************
 fit.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <math.h>
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Constants.H"
#include "BOOM/GSL/Optimizer.H"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
using namespace std;
using namespace BOOM;

const double EPSILON=0.1; 
const double stepSize=0.001;
const double tolerance=0.1;
const double gradientThreshold=0.001;
const GSL::StoppingCriterion criterion=GSL::BY_EITHER;

typedef pair<double,double> Point;
typedef Vector<Point> Curve;


/****************************************************************
                        enum DistributionType
 ****************************************************************/
enum DistributionType {
  DISTR_LOGNORMAL,
  DISTR_GUMBEL,
  DISTR_PARETO,
  DISTR_BETA,
  DISTR_STUDENT,
  DISTR_F,
  DISTR_CHI,
  DISTR_GAMMA,
  DISTR_RAYLEIGH,
  DISTR_WEIBULL,
  DISTR_EXPONENTIAL,
  DISTR_LAPLACE,
  DISTR_CAUCHY,
  DISTR_LOGISTIC
};
ostream &operator<<(ostream &,DistributionType);
DistributionType distributionFromString(const String &);
int numParms(DistributionType);


/****************************************************************
                          enum ParmIndex
 ****************************************************************/
enum ParmIndex {
  PARM_ORIGIN, // always used
  PARM_SCALE,  // always used
  PARM_LAMBDA, // always used
  PARM_MU,     // used by all 1-parm distributions
  PARM_SIGMA,  // used by all 2-parm distributions
};
const int MAX_PARM_INDEX=4;
ostream &operator<<(ostream &,ParmIndex);
inline ParmIndex &operator++(ParmIndex &i) {i=(ParmIndex)i+1;}


/****************************************************************
                       class ObjectiveFunc
 ****************************************************************/
struct ObjectiveFunc : public GSL::ObjectiveFunction {
  ObjectiveFunc(Curve &targetCurve,DistributionType,int DIM,
		int NUM_ITERATIONS,int MAX_STEPS,int CURVE_LEN);
  virtual double f(const GSL::Vector &currentPoint);
  virtual void gradient(const GSL::Vector &currentPoint,
			GSL::Vector &gradient);
  int DIM;
  int NUM_ITERATIONS;
  int MAX_STEPS;
  int CURVE_LEN;
  DistributionType distrType;
  Curve &targetCurve;
  void unpack(Vector<double> &parms,double &mu,double &sigma,
	      double &lambda,double &scale,double &origin);
  double evaluatePoint(Vector<double> &parms,Curve &);
  void buildCurve(Vector<double> &parms,Curve &,int len);
  double diff(Curve &points1,Curve &points2);
  void dumpParms(Vector<double> &parms,ostream &);
};


/****************************************************************
                        class Application
 ****************************************************************/
class Application {
  ObjectiveFunc *objectiveFunc;
  bool HILL_CLIMB;
  void optimize(Curve &,Vector<double> &bestParms,double &bestVal);
  void hillClimber(Curve &c,Vector<double> &bestParms,
		   double &bestVal,Vector<double> deltas);
public:
  Application();
  int main(int argc,char *argv[]);
  static GSL::Vector *toGSL(const Vector<double> &);
  static Vector<double> *fromGSL(const GSL::Vector &);
};


/****************************************************************
                              main()
 ****************************************************************/
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




/****************************************************************
                       Application methods
 ****************************************************************/

Application::Application()
  : HILL_CLIMB(false)
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"i:s:");
  if(cmd.numArgs()!=4)
    throw String("\n\
Fits lambda*P(length>x) to a given curve\n\
\n\
usage:  fit [options] <in.hist> <out.hist> <out.parms> <distribution>\n\
  where distribution is:\n\
     lognormal|gumbel|pareto|beta|student|F|chi|gamma|\n\
     rayleigh|weibull|exponential|laplace|cauchy|logistic\n\
  -i : max iterations\n\
  -s : max steps in line search\n\
");
  String infile=cmd.arg(0);
  String outfile=cmd.arg(1);
  String parmFile=cmd.arg(2);
  String distrName=cmd.arg(3);
  DistributionType distrType=distributionFromString(distrName);
  int DIM=3+numParms(distrType);
  int NUM_ITERATIONS=10, MAX_STEPS=100;
  if(cmd.option('i')) NUM_ITERATIONS=cmd.optParm('i').asInt();
  if(cmd.option('s')) MAX_STEPS=cmd.optParm('s').asInt();

  // Load the target curve
  Curve targetCurve;
  ifstream is(infile.c_str());
  while(!is.eof()) {
    double x, y;
    is>>x;
    if(is.eof()) break;
    if(x<1) continue;
    is>>y;
    targetCurve.push_back(Point(x,y));
  }
  is.close();
  int CURVE_LEN=targetCurve.size();

  // Set up initial point
  Vector<double> bestParms(MAX_PARM_INDEX+1);
  bestParms[PARM_MU]=0.1;
  bestParms[PARM_SIGMA]=0.9;
  bestParms[PARM_LAMBDA]=0.5;
  bestParms[PARM_SCALE]=400;
  bestParms[PARM_ORIGIN]=0;

  // Set up objective function
  objectiveFunc=new ObjectiveFunc(targetCurve,distrType,DIM,NUM_ITERATIONS,
				  MAX_STEPS,CURVE_LEN);

  // Perform optimization
  double bestVal;
  optimize(targetCurve,bestParms,bestVal);

  // Emit results
  cout<<"BEST: ";objectiveFunc->dumpParms(bestParms,cout);
  cout<<" => "<<bestVal<<endl;
  ofstream osParms(parmFile.c_str());
  for(ParmIndex i=0 ; i<=MAX_PARM_INDEX ; ++i) 
    osParms<<i<<"="<<bestParms[i]<<endl;
  osParms.close();
  Curve c;
  objectiveFunc->buildCurve(bestParms,c,CURVE_LEN);
  int n=c.size();
  ofstream os(outfile.c_str());
  for(int i=0 ; i<n ; ++i) {
    int x=i+1;
    double y=c[i].second;
    os<<x<<"\t"<<y<<endl;
  }
  return 0;
}



void Application::optimize(Curve &c,Vector<double> &bestParms,
			   double &bestVal)
{
  if(HILL_CLIMB) {
    Vector<double> deltas(MAX_PARM_INDEX+1);
    deltas[PARM_MU]=0.01174; 
    deltas[PARM_SIGMA]=0.01234;
    deltas[PARM_LAMBDA]=0.0126371;
    deltas[PARM_SCALE]=10.82761;
    deltas[PARM_ORIGIN]=10.17621;
    hillClimber(c,bestParms,bestVal,deltas);
  }
  else {
    GSL::Vector *initialPoint=toGSL(bestParms);
    cout<<initialPoint->getDim()<<" "<<stepSize<<" "<<tolerance<<" "<<int(criterion)<<" "<<objectiveFunc->NUM_ITERATIONS<<" "<<int(GSL::BFGS)<<endl;
    GSL::Optimizer O(GSL::BFGS,*objectiveFunc,*initialPoint,stepSize,criterion,
		     tolerance,gradientThreshold,
		     objectiveFunc->NUM_ITERATIONS);
    bool status=O.run();
    if(!status) cout<<"Error during optimization"<<endl;
    const GSL::Vector &opt=O.getOptimalPoint();
    Vector<double> *tmp=fromGSL(opt);
    bestParms=*tmp;
    delete tmp;
    delete initialPoint;
  }
}



GSL::Vector *Application::toGSL(const Vector<double> &v)
{
  int n=v.size();
  GSL::Vector &g=*new GSL::Vector(n);
  for(int i=0 ; i<n ; ++i) g[i]=v[i];
  return &g;
}



Vector<double> *Application::fromGSL(const GSL::Vector &v)
{
  int n=v.getDim();
  Vector<double> &g=*new Vector<double>(n);
  for(int i=0 ; i<n ; ++i) g[i]=v[i];
  return &g;
}



void Application::hillClimber(Curve &c,Vector<double> &bestParms,
				double &bestVal,Vector<double> deltas)
{
  Vector<double> currentPt=bestParms;
  double currentVal=objectiveFunc->evaluatePoint(currentPt,c);
  cout<<"initial score: ";
  objectiveFunc->dumpParms(currentPt,cout);
  cout<<" => "<<currentVal<<endl;
  for(int iter=0 ; iter<objectiveFunc->NUM_ITERATIONS ; ++iter) {
    cout<<"iteration #"<<iter<<endl;
    // Iterate through all parameters (dimensions)
    for(int dim=objectiveFunc->DIM-1 ; dim>=0 ; --dim) {
      int steps=1;
      cout<<"optimizing dimension #"<<dim<<endl;
      // Perform a line search in this dimension
      int sign=1;
      double delta=deltas[dim]*sign;
      while(1) { // While accuracy continues to improve...
	if(++steps>objectiveFunc->MAX_STEPS) break;

	// Take a step along one dimension
	Vector<double> newPt=currentPt;
	newPt[dim]+=delta;
	
	// Evaluate accuracy at new point
	double newVal=objectiveFunc->evaluatePoint(newPt,c);
	cout.precision(12);
	//for(int j=0 ; j<newPt.size() ; ++j) cout<<"\t"<<newPt[j];
	objectiveFunc->dumpParms(newPt,cout);
	cout<<" => "<<newVal;
        
	// Decide whether to continue or turn back
	if(newVal>currentVal)
	  { currentPt=newPt; currentVal=newVal; cout<<"\tACCEPT"<<endl;}
	else {
	  cout<<"\tREJECT"<<endl;
	  if(sign==-1) break;
	  sign*=-1;
	  delta=deltas[dim]*sign;
	}
      }
    }
  }
  bestVal=currentVal;
  bestParms=currentPt;
}



/****************************************************************
                      ObjectiveFunc methods
 ****************************************************************/

ObjectiveFunc::ObjectiveFunc(Curve &targetCurve,DistributionType d,
			     int DIM,int NUM_ITERATIONS,int MAX_STEPS,
			     int CURVE_LEN)
  : targetCurve(targetCurve), distrType(d), DIM(DIM), 
    NUM_ITERATIONS(NUM_ITERATIONS), MAX_STEPS(MAX_STEPS), 
    CURVE_LEN(CURVE_LEN)
{
  // ctor
}



void ObjectiveFunc::unpack(Vector<double> &parms,double &mu,double &sigma,
			 double &lambda,double &scale,double &origin)
{
  mu=parms[PARM_MU];
  sigma=parms[PARM_SIGMA];
  lambda=parms[PARM_LAMBDA];
  scale=parms[PARM_SCALE];
  origin=parms[PARM_ORIGIN];
}


double ObjectiveFunc::evaluatePoint(Vector<double> &parms,Curve &c)
{
  double mu, sigma, lambda, scale, origin;
  unpack(parms,mu,sigma,lambda,scale,origin);
  if(lambda<0 || mu<0 || scale<0) return SMALLEST_DOUBLE;
  Curve newCurve;
  buildCurve(parms,newCurve,CURVE_LEN);
  return diff(c,newCurve);
}



void ObjectiveFunc::buildCurve(Vector<double> &parms,Curve &curve,int L)
{
  double mu, sigma, lambda, scale, origin;
  unpack(parms,mu,sigma,lambda,scale,origin);
  for(int x=1 ; x<=L ; ++x) {
    double dist=(x-origin)/scale;
    double p;
    switch(distrType) {
    case DISTR_LOGNORMAL:      
      p=lambda*gsl_cdf_lognormal_Q(dist,mu,sigma);
      //double sig2=sigma*sigma;
      //p=lambda*(.5-.5*erf((log(dist)-mu)/sqrt(2*sig2)));
      break;
    case DISTR_GUMBEL:     
      p=lambda*gsl_cdf_gumbel2_Q(dist,mu,sigma);
      break;
    case DISTR_PARETO:     
      p=lambda*gsl_cdf_pareto_Q(dist,mu,sigma);
      //p=lambda*(dist>=sigma ? pow(sigma/dist,mu) : 1.0);
      break;
    case DISTR_BETA:       
      p=lambda*gsl_cdf_beta_Q(dist,mu,sigma);
      break;
    case DISTR_STUDENT:   
      p=lambda*gsl_cdf_tdist_Q(dist,mu);
      break;
    case DISTR_F:          
      p=lambda*gsl_cdf_fdist_Q(dist,mu,sigma);
      break;
    case DISTR_CHI:        
      p=lambda*gsl_cdf_chisq_Q(dist,mu);
      break;
    case DISTR_GAMMA:      
      p=lambda*gsl_cdf_gamma_Q(dist,mu,sigma);
      break;
    case DISTR_RAYLEIGH:   
      p=lambda*gsl_cdf_rayleigh_Q(dist,mu);
      break;
    case DISTR_WEIBULL:    
      p=lambda*gsl_cdf_weibull_Q(dist,mu,sigma);
      break;
    case DISTR_EXPONENTIAL:
      p=lambda*gsl_cdf_exponential_Q(dist,mu);
      break;
    case DISTR_LAPLACE:    
      p=lambda*gsl_cdf_laplace_Q(dist,mu);
      break;
    case DISTR_CAUCHY:     
      p=lambda*gsl_cdf_cauchy_Q(dist,mu);
      break;
    case DISTR_LOGISTIC:   
      p=lambda*gsl_cdf_logistic_Q(dist,mu);
      break;
    default: INTERNAL_ERROR;
    }
    curve.push_back(Point(x,p));
  }
}



double ObjectiveFunc::diff(Curve &points2,Curve &points1) {
  int n1=points1.size(), n2=points2.size();
  int n=max(n1,n2);
  double d=0.0;
  for(int i=0 ; i<n; ++i) {
    double y1=i<n1 ? points1[i].second : 0.0;
    double y2=i<n2 ? points2[i].second : 0.0;
    d+=fabs(y1-y2);
  }
  return d;
}



ostream &operator<<(ostream &os,DistributionType d)
{
  switch(d) {
  case DISTR_LOGNORMAL:  os<<"lognormal"<<endl; break;
  case DISTR_GUMBEL:     os<<"gumbel"<<endl; break;
  case DISTR_PARETO:     os<<"pareto"<<endl; break;
  case DISTR_BETA:       os<<"beta"<<endl; break;
  case DISTR_STUDENT:    os<<"student"<<endl; break;
  case DISTR_F:          os<<"F"<<endl; break;
  case DISTR_CHI:        os<<"chi"<<endl; break;
  case DISTR_GAMMA:      os<<"gamma"<<endl; break;
  case DISTR_RAYLEIGH:   os<<"rayleigh"<<endl; break;
  case DISTR_WEIBULL:    os<<"weibull"<<endl; break;
  case DISTR_EXPONENTIAL:os<<"exponential"<<endl; break;
  case DISTR_LAPLACE:    os<<"laplace"<<endl; break;
  case DISTR_CAUCHY:     os<<"cauchy"<<endl; break;
  case DISTR_LOGISTIC:   os<<"logistic"<<endl; break;
  default: INTERNAL_ERROR;
  }
};



DistributionType distributionFromString(const String &s)
{
  if(s=="lognormal") return DISTR_LOGNORMAL;
  if(s=="gumbel") return DISTR_GUMBEL;
  if(s=="pareto") return DISTR_PARETO;
  if(s=="beta") return DISTR_BETA;
  if(s=="student") return DISTR_STUDENT;
  if(s=="F") return DISTR_F;
  if(s=="chi") return DISTR_CHI;
  if(s=="gamma") return DISTR_GAMMA;
  if(s=="rayleigh") return DISTR_RAYLEIGH;
  if(s=="weibull") return DISTR_WEIBULL;
  if(s=="exponential") return DISTR_EXPONENTIAL;
  if(s=="laplace") return DISTR_LAPLACE;
  if(s=="cauchy") return DISTR_CAUCHY;
  if(s=="logistic") return DISTR_LOGISTIC;
  throw s+": unknown distribution";
}



ostream &operator<<(ostream &os,ParmIndex p) {
  switch(p) {
  case PARM_MU:     os<<"mu"; break;
  case PARM_SIGMA:  os<<"sigma"; break;
  case PARM_LAMBDA: os<<"lambda"; break;
  case PARM_SCALE:  os<<"scale"; break;
  case PARM_ORIGIN: os<<"origin"; break;
  default: INTERNAL_ERROR;
  }
};



int numParms(DistributionType d)
{
  switch(d) {
  case DISTR_LOGNORMAL: return 2;
  case DISTR_GUMBEL: return 2;
  case DISTR_PARETO: return 2;
  case DISTR_BETA: return 2;
  case DISTR_STUDENT: return 1;
  case DISTR_F: return 2;
  case DISTR_CHI: return 1;
  case DISTR_GAMMA: return 2;
  case DISTR_RAYLEIGH: return 1;
  case DISTR_WEIBULL: return 2;
  case DISTR_EXPONENTIAL: return 1;
  case DISTR_LAPLACE: return 1;
  case DISTR_CAUCHY: return 1;
  case DISTR_LOGISTIC: return 1;
  default: INTERNAL_ERROR;
  }
}



void ObjectiveFunc::dumpParms(Vector<double> &parms,ostream &os)
{
  os.precision(12);
  for(ParmIndex i=0 ; i<DIM ; ++i) {
    os<<i<<"="<<parms[i];
    if(i<DIM) os<<",";
  }
}



double ObjectiveFunc::f(const GSL::Vector &currentPoint)
{
  Vector<double> *parms=Application::fromGSL(currentPoint);
  double v=evaluatePoint(*parms,targetCurve);
  delete parms;
  cout<<"objective: "<<v<<endl;
  return v;
}



void ObjectiveFunc::gradient(const GSL::Vector &currentPoint,
			     GSL::Vector &gradient)
{
  const double epsilon=EPSILON;
  GSL::Vector perturbed=currentPoint;
  int n=gradient.getDim();
  for(int i=0 ; i<n ; ++i)
    {
      double &x=perturbed[i];
      const double x0=x;
      double dx=epsilon*x;
      if(dx<epsilon) dx=epsilon;
      double temp=x+dx; 
      dx=temp-x; // make sure it's an "exact value" in the machine
      double twoDX=2*dx;
      x-=dx;
      double y1=f(perturbed);
      x+=twoDX;
      double y2=f(perturbed);
      gradient[i]=(y2-y1)/twoDX;
      x=x0;
      //if(!isFinite(gradient[i])) cout<<"x0="<<x0<<" dx="<<dx<<" twoDX="<<twoDX<<" y1="<<y1<<" y2="<<y2<<endl;
    }
  cout<<"gradient = "<<gradient<<endl;
}







