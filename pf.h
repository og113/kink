//parameters and functions for pi.cc
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <ctype.h>
#include <limits>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "gnuplot_i.hpp"

using namespace std;

typedef unsigned int uint;
typedef unsigned long int lint;
typedef complex<double> comp;
typedef vector<unsigned int> intVec;
typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;

complex<double> ii(0.0,1.0);
#define pi 3.14159265359

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//primary parameters
unsigned int N;
unsigned int Na;
unsigned int Nb;
unsigned int Nc;
double LoR; //L/R
double dE;
double Tb;  //b section includes both corner points
double theta;

//derived parameters
unsigned int NT;
double epsilon;
double R; //size of bubble
double Gamma; //equals exp(-theta)
double angle; //not a primary parameter, just used to make L
double L;
double a; //step sizes in each spatial dimension
double b; //step sizes in time
double Ta;
double Tc;
vector<double> minima(2);
double mass2; //as derived from V''
double A; //gives parameter in V2, equal to 0.4 in DL

//determining number of runs
double closenessA; //action
double closenessS; //solution (i.e. minusDS)
double closenessSM; //solution max
double closenessD; //delta
double closenessC; //calculation
double closenessCon; //energy conservation
double closenessL; //linearisation of energy
double closenessT; //true energy versus linear energy
double closenessP; //checking lattice small enough for momenta
double closenessR; //regularization
double closenessIE; //imaginary part of the energy
double closenessCL; //continuum approx to linear erg and num
double closenessON; //on shell vs off shell linear energy
double closenessAB; //a_k=Gamma*b*_k
double closenessLR; //p(i,j)->p_lin(i,j) as i->0
double closenessABNE; //linNum and linErg as calculated by a_k*b*_k agree with other calcs

//parameters determining input phi
//struct to hold answers to questions
struct aqStruct
	{
	string inputChoice;
	string inputTimeNumber;
	string inputLoop;
	double maxTheta;
	unsigned int totalLoops;
	string loopChoice;
	double minValue;
	double maxValue;
	string printChoice;
	unsigned int printRun;
	};
aqStruct aq; //struct to hold user responses
double reg; //small parameter multiplying regulatory term
unsigned int cutoff; //cutoff for momentum integrations
string inP; //b for bubble, p for periodic instaton, f for from file
string pot; //pot[0] gives 1 or 2, pot[1] gives r (regularised) or n (not)
string inF; //file input from where, m for main, p for pi
int firstLoop, lastLoop; //first and last loop to load from, not unsigned for comparison later
double alpha; //gives span over which tanh is used
double open; //value of 0 assigns all weight to boundary, value of 1 to neighbour of boundary
double amp; //ammount of negative eigenvector added to bubble for Tb>R
double negVal; //the negative eigenvalue
unsigned int negEigDone; //has the negEig been found before? 1 if yes, 0 if no
string zmt; //dictates how time zero mode is dealt with
string zmx; //dictates how x zero mode is dealt with
string bds; //dictates how boundaries are dealt with in main
double epsilon0; //the value of epsilon when the minima are degenerate
double S1;
double r0; //the minimum radius if pot[0]=='3'

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//simple generic functions

//factorial function
int factorial(const int& f_input)
	{
	int f_result = 1;
	for (int l=0; l<f_input; l++) f_result *= (l+1);
	return f_result;
	}
	
//gives absolute value of a number
double absolute (const double& amplitude)
	{
	double abs_amplitude;
	if (amplitude > 0)
		{
		abs_amplitude = amplitude;
		}
	else
		{
		abs_amplitude = -amplitude;
		}
	return abs_amplitude;
	}
	
//gives absolute measure of difference between two numbers
double absDiff (const double& numA, const double& numB) {
	return 2.0*abs(numA-numB)/sqrt(numA*numA+numB*numB);
}
double absDiff (const comp& numA, const comp& numB) {
	return 2.0*abs(numA-numB)/sqrt(norm(numA)+norm(numB));
}

double absDiff(const vec& vecA, const vec& vecB) {
	double normSqrdA = vecA.squaredNorm();
	double normSqrdB = vecB.squaredNorm();
	vec diff = vecA-vecB;
	double normDiff = diff.norm();
	return 2.0*normDiff/pow(normSqrdA+normSqrdB,0.5);
}

double absDiff(const cVec& vecA, const cVec& vecB) {
	double normSqrdA = vecA.squaredNorm();
	double normSqrdB = vecB.squaredNorm();
	cVec diff = vecA-vecB;
	double normDiff = diff.norm();
	return 2.0*normDiff/pow(normSqrdA+normSqrdB,0.5);
}
	
//function giving location of smallest element of a vector of type T
template <typename T>
unsigned int smallestLoc(const vector <T> & inVector)
	{
	unsigned int loc = 0;
	for(unsigned int l=1;l<inVector.size();l++)
		{
		if (absolute(inVector[l])<absolute(inVector[loc]))
			{
			loc = l;
			}
		}
	return loc;
	}
	
//count non-empty lines of a file
unsigned int countLines(const string & file_to_count)
	{
	ifstream fin;
	fin.open(file_to_count.c_str());
	if (!fin.good()) cerr << "countLines error: " << file_to_count << " not opened properly." << endl;
	string line;
	unsigned int counter = 0;
	while(!fin.eof())
		{
		getline(fin,line);
		if(line.empty()) continue;
		counter++;
		}		
	fin.close();
    return counter;
	}
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//generic gsl derived functions

//function to find a root of function FDF, given initial guess, using newton method
double rootFinder(gsl_function_fdf * xFDF, double rootGuess)
	{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	double x0, x = rootGuess;

	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc (T);
	gsl_root_fdfsolver_set (s, xFDF, x);

	do
		{
		  iter++;
		  status = gsl_root_fdfsolver_iterate (s);
		  x0 = x;
		  x = gsl_root_fdfsolver_root (s);
		  status = gsl_root_test_delta (x, x0, 0, DBL_MIN);
		}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fdfsolver_free (s);
	return x;
	}
	
//function to find a root of function FDF, given initial guess, and lower and upper bounds, using brent method
double brentRootFinder(gsl_function * xF, const double & rootGuess, const double & rootLower, const double &rootUpper)
	{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double x = rootGuess;
	double x_lo = rootLower;
	double x_hi = rootUpper;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, xF, x_lo, x_hi);

	do
		{
		  iter++;
		  status = gsl_root_fsolver_iterate (s);
		  x = gsl_root_fsolver_root (s);
		  x_lo = gsl_root_fsolver_x_lower (s);
		  x_hi = gsl_root_fsolver_x_upper (s);
      	  status = gsl_root_test_interval (x_lo, x_hi,
                                       0, DBL_MIN);
		}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
	return x;
	}
	
//function to find the minimum of a gsl function F
double brentMinimum (gsl_function * xF, const double & minimumGuess, const double & minimumLower, const double & minimumUpper)
	{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m = minimumGuess;
	double tempLower, tempUpper;

	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc (T);
	gsl_min_fminimizer_set (s, xF, m, minimumLower, minimumUpper);
	
	do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      tempLower = gsl_min_fminimizer_x_lower (s);
      tempUpper = gsl_min_fminimizer_x_upper (s);

      status = gsl_min_test_interval (tempLower, tempUpper, 1.0e-12, DBL_MIN);

      if (status == GSL_SUCCESS)
      	{
        //printf ("brentMinimum converged\n");
        }
    }
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_min_fminimizer_free (s);
	
	return m;
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//potential functions and their derivatives

//structs containing parameters
struct params_for_V {double epsi; double aa;};
struct params_for_V paramsV, paramsV0;
struct void_params {};
struct void_params paramsVoid;

//////////////////////////////potentials//////////////////////////////
//V1
//NB - gsl doesn't like const T &; default parameters set
template <class T> T V1 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	return pow(pow(phi,2)-1.0,2.0)/8.0 - epsi*(phi-1.0)/2.0;
	}

comp V1c (const comp phi) { return V1(phi); }

//Z, for V2
template <class T> T Z (const T phi)
	{
	return exp(-pow(phi,2.0))*(phi + pow(phi,3.0) + pow(phi,5.0));
	}
	
//V2
template <class T> T V2 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	double aa = (params->aa);
	return 0.5*pow(phi+1.0,2.0)*(1.0-epsi*Z((phi-1.0)/aa));
	}
comp V2c (const comp phi) { return V2(phi); }

//V3
template <class T> T V3 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi); //epsi here takes the role of r
	return pow(phi,2.0)/2.0 - pow(phi,4.0)/4.0/pow(epsi,2.0);
	}
comp V3c (const comp phi) { return V3(phi); }


//Vr	
comp VrFn (const comp & phi, const double & minimaL, const double & minimaR)
	{
	return pow(phi-minimaL,4.0)*pow(phi-minimaR,4.0)/4.0;
	}

//declaration of function pointer V
comp (*V) (const comp phi); //a function pointer
double (*Vd) (const double phi, void * parameters);

//////////////////////////////first derivatives of potentials//////////////////////////////
//dV1
template <class T> T dV1 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	return phi*(pow(phi,2)-1.0)/2.0 - epsi/2.0;
	}	
comp dV1c (const comp phi) { return dV1(phi); }

//dZ for dV2
template <class T> T dZ (const T phi)
	{
	return exp(-pow(phi,2.0))*( 1.0 + pow(phi,2.0) + 3.0*pow(phi,4.0) - 2.0*pow(phi,6.0) );
	}
	
//dV2
template <class T> T dV2 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	double aa = (params->aa);
	return (phi+1.0)*(1.0-epsi*Z((phi-1.0)/aa)) - 0.5*pow(phi+1.0,2.0)*(epsi/aa)*dZ((phi-1.0)/aa);
	}
comp dV2c (const comp phi) { return dV2(phi); }

//dV3
template <class T> T dV3 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	return phi - pow(phi,3.0)/pow(epsi,2.0);
	}	
comp dV3c (const comp phi) { return dV3(phi); }
	
//dVr
comp dVrFn (const comp & phi, const double & minimaL, const double & minimaR)
	{
	return pow(phi-minimaL,3.0)*pow(phi-minimaR,4.0) + pow(phi-minimaL,4.0)*pow(phi-minimaR,3.0);
	}
	
comp (*dV) (const comp phi);
double (*dVd) (const double phi, void * parameters);
	
//////////////////////////////second derivatives of potentials//////////////////////////////
//ddV1
template <class T> T ddV1 (const T phi, void * parameters = &paramsVoid)
	{
	return (3.0*pow(phi,2)-1.0)/2.0;
	}
comp ddV1c (const comp phi) { return ddV1(phi); }

//ddZ for ddV2
template <class T> T ddZ (const T phi)
	{
	return exp(-pow(phi,2.0))*2.0*pow(phi,3.0)*(5.0 - 9.0*pow(phi,2.0) + 2.0*pow(phi,4.0));
	}

//ddV2
template <class T> T ddV2 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	double aa = (params->aa);
	return (1.0-epsi*Z((phi-1.0)/aa)) - (phi+1.0)*(epsi/aa)*dZ((phi-1.0)/aa)\
					+ 0.5*pow(phi+1.0,2.0)*(epsi/pow(aa,2.0))*ddZ((phi-1.0)/aa);
	}
comp ddV2c (const comp phi) { return ddV2(phi); }

//ddV3
template <class T> T ddV3 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	return 1.0 - 3.0*pow(phi,2.0)/pow(epsi,2.0);
	}
comp ddV3c (const comp phi) { return ddV3(phi); }

//ddVr	
comp ddVrFn (const comp & phi, const double & minimaL, const double & minimaR)
	{
	return 3.0*pow(phi-minimaL,2.0)*pow(phi-minimaR,4.0) + 8.0*pow(phi-minimaL,3.0)*pow(phi-minimaR,3.0)\
				+ 3.0*pow(phi-minimaL,4.0)*pow(phi-minimaR,2.0);
	}

comp (*ddV) (const comp phi);
double (*ddVd) (const double phi, void * parameters);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//extra functions
	
//V and dV FDF gsl functions	
void VdV (double x, void * parameters, double * f, double* df) 
	{
	*f =  Vd(x,parameters);
	*df = dVd(x,parameters);
	}
	
void dVddV (double x, void * parameters, double * f, double* df) 
	{
	*f =  dVd(x,parameters);
	*df = ddVd(x,parameters);
	}
	
//energy change gsl function : V(minima[1])-V(minima[0])-dE
struct ec_params {double aa; double minima0; double minima1; double de; };

double ec (double epsi, void * parameters)
	{
	struct ec_params * paramsIn = (struct ec_params *)parameters;
	struct params_for_V paramsOut;
	paramsOut.epsi = epsi;
	paramsOut.aa = (paramsIn->aa);
	double minima1 = (paramsIn->minima1);
	double minima0 = (paramsIn->minima0);
	double de = (paramsIn->de);
	return Vd(minima0,&paramsOut) - Vd(minima1,&paramsOut) - de;
	}
	
//S1 integrand
double s1Integrand (double x, void * parameters)
	{
	return pow(2.0*Vd(x,parameters),0.5);
	}

//rho integrand
double rhoIntegrand (double x, void * parameters)
	{
	return pow(2.0*Vd(x,parameters),-0.5);
	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//functions to calculate roots and minima
	
//function to give the three roots of FDF given lower and upper limits on them and a number of loops to try
vector <double> rootsFn (gsl_function_fdf * xFDF, const double & lowLimit, const double & highLimit,\
							const unsigned int & rootLoops)
	{
	vector <double> roots;
	for (unsigned int j=0;j<rootLoops;j++)
		{
		double x = lowLimit+(absolute(highLimit)+absolute(lowLimit))*j/(rootLoops-1.0);
		x = rootFinder(xFDF,x);
	
		if (j==0)
			{
			roots.push_back(x);
			}
		else
			{
			unsigned int test = 0;
			for (unsigned int k=0;k<roots.size();k++)
				{ 
				if (absolute(x-roots[k])>1.0e-6)
					{
					test++;
					}
				}
			if (test==roots.size())
				{
				roots.push_back(x);
				}
			}
		}
	
		if (roots.size()!=3)
			{
			cout << "rootsFn error: only found " << roots.size() << " roots, not 3" << endl;
			}
	return roots;
	}
	
//program to find epsilon given gsls function df and dE
void epsilonFn (gsl_function * xF, gsl_function * xEC, double * xdE, double * xEpsilon, vector<double>* xMinima)
	{
	double closenessdE = 1.0e-14;
	vector<double> dE_test(1);	dE_test[0] = 1.0;
	double newdE = dE;
	struct params_for_V * Fparameters = (struct params_for_V *) (*xF).params;
	struct ec_params * ECparameters = (struct ec_params *) (*xEC).params;
	unsigned int counter = 0;
	unsigned int maxCounter = 1e4;
	while (dE_test.back()>closenessdE)
		{
		//find roots of ec(epsilon)=0
		*xEpsilon = brentRootFinder(xEC,*xEpsilon,*xEpsilon/2.0,*xEpsilon*2.0);
		//assign new value of epsilon to xF
		(*Fparameters).epsi = *xEpsilon;
		(*xF).params = Fparameters;
		//finding new roots of dV(phi)=0
		(*xMinima)[0] = brentMinimum(xF,-1.0,-3.0,0.0);
		(*xMinima)[1] = brentMinimum(xF,1.2,0.5,3.0);
		//assign new roots to xECDF
		(*ECparameters).minima0 = (*xMinima)[0];
		(*ECparameters).minima1 = (*xMinima)[1];
		(*xEC).params = ECparameters;
		//evaluating new dE
		newdE = (*(*xEC).function)(*xEpsilon,ECparameters) + *xdE;
		//evaluating test
		if (absolute(*xdE)>1.0e-16)
			{
			dE_test.push_back(absolute((newdE-(*xdE))/(*xdE)));
			}
		else
			{
			dE_test.push_back(absolute(newdE-(*xdE)));
			}
		counter++;
		//test if too many runs
		if (counter>maxCounter)
			{
			cout << "epsilonFn error, more that " << maxCounter << " loops, consider reducing closenessdE" << endl;
			cout << "dE_test.back() = " << dE_test.back() << " , closenessdE = " << closenessdE << endl;
			cout << "dE = " << *xdE << " , minima[0] = " << (*xMinima)[0] << " , minima[1] = " << (*xMinima)[1];
			cout << " , epsilon = " << *xEpsilon << endl << endl;
			break;
			}
		}
	*xdE = newdE;
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//lattice functions

//inverse of intCoord
unsigned int c(const unsigned int& t, const unsigned int& x, const unsigned int& Nt) {
	return t + x*Nt;
}

//function which find integer coordinates in 2d
unsigned int intCoord(const unsigned int& locNum, const int& direction, const unsigned int& xNt)
	{
	unsigned int XintCoord;
	unsigned int x = floor(locNum/xNt);
	if (direction==1) 	XintCoord = x;
	else 				XintCoord = locNum - x*xNt;
	return XintCoord;	
	}
	
//simple time
comp simpleTime (const unsigned int& time, const unsigned int& the_Na=Na, const unsigned int& the_Nb=Nb, const unsigned int& the_Nc=Nc, const double& the_Tb = Tb, const double& the_b=b)
	{
	comp xTime;
	if ( time < the_Na)
		{
		double temp = (double)time;
		temp -= (double)the_Na;
		xTime = the_b*temp + ii*the_Tb;
		}
	else if (time < (the_Na+the_Nb))
		{
		double temp = (double)time;
		temp -= the_Na; //as complex doesn't support (complex double)*integer (though it does support double*integer added to a complex double) - and as int to double seems to cock up here (perhaps because the integers are unsigned)
		xTime = ii*(the_Tb - the_b*temp);
		}
	else
		{
		double temp = (double)time;
		temp -= (double)the_Na;
		temp -= (double)the_Nb;
		xTime = the_b*(temp+1.0); //the 1.0 is because the corner is part of the vertical contour
		}
	return xTime;
	}
	
//simple space functions
double (*simpleSpace) (const unsigned int& space, const double& the_L, const double& the_a);	

//simple space for a box
double simpleSpaceBox (const unsigned int& space, const double& the_L, const double& the_a)
	{
	return -L/2.0 + space*the_a;;
	}
	
//simple space for a sphere
double simpleSpaceSphere (const unsigned int& space, const double& the_L=L, const double& the_a=a)
	{
	return r0+space*the_a;
	}
	
//gives values of coordinates in whole spacetime
comp coord(const unsigned int& locNum,const int& direction, const unsigned int& Nt=NT)
	{
	if (direction==0)		return simpleTime(intCoord(locNum,0,NT));
	else if (direction==1)	return simpleSpace(intCoord(locNum,1,NT),L,a);
	else return 0.0;
	}

//gives values of coordinates on section AB
complex<double> coordA(const unsigned int& locNum,const int& direction)
	{
	complex<double> XcoordA;
	if (direction==0)
		{
		unsigned int t = intCoord(locNum,0,Na);
		double temp = (double)t;
		temp -= (double)Na;
		XcoordA = b*temp + ii*Tb;
		}
	if (direction==1)
		{
		unsigned int x = intCoord(locNum,1,Na);
		XcoordA = simpleSpace(x,L,a);
		}
	return XcoordA;
	}

//gives values of coordinates on section BC
complex<double> coordB(const unsigned int& locNum,const int& direction)
	{
	complex<double> XcoordB;
	if (direction==0)
		{
		unsigned int t = intCoord(locNum,0,Nb);
		double temp = (double)t;
		XcoordB = ii*(Tb - b*temp);
		}
	if (direction==1)
		{
		unsigned int x = intCoord(locNum,1,Nb);
		XcoordB = simpleSpace(x,L,a);
		}
	return XcoordB;
	}
	
//gives values of coordinates on section CD
complex<double> coordC(const unsigned int& locNum,const int& direction)
	{
	complex<double> XcoordC;
	if (direction==0)
		{
		unsigned int t = intCoord(locNum,0,Nc);
		XcoordC = b*t;
		}
	if (direction==1)
		{
		unsigned int x = intCoord(locNum,1,Nc);
		XcoordC = simpleSpace(x,L,a);
		}
	return XcoordC;
	}

long int periodic(const lint& locNum, const unsigned int& direction, const signed int& sign, const unsigned int& xNt, const unsigned int & xNx) //periodic in space but not time, degree refers to the number of neighbours, 1 is for just positive neighbours, 2 is for both
	{
	long int neighLocation = -1; //this is the result if there are no neighbours for the given values of the argument
	unsigned int c = intCoord(locNum,direction,xNt);
	if (direction==0)
		{
		if (sign==1 and c!=(xNt-1))
			{
			neighLocation = locNum+1;
			}
		else if (sign==-1 and c!=0)
			{
			neighLocation = locNum-1;
			}
		}
	else if (c==0 and sign==-1)
		{
		neighLocation = locNum+(xNx-1)*xNt;
		}
	else if (c==(xNx-1) and sign==1)
		{
		neighLocation = locNum-(xNx-1)*(int)xNt;
		}
	else
		{
		neighLocation = locNum+sign*(int)xNt;
		}
	return neighLocation;
	}
	
//spherical system, reflective at r=0, nothing at r=R
long int spherical(const lint& locNum, const unsigned int& direction, const signed int& sign, const unsigned int& xNt, const unsigned int & xNx) //periodic in space but not time, degree refers to the number of neighbours, 1 is for just positive neighbours, 2 is for both
	{
	long int neighLocation = -1; 
	unsigned int c = intCoord(locNum,direction,xNt);
	if (direction==0)
		{
		if (sign==1 and c!=(xNt-1))
			{
			neighLocation = locNum+1;
			}
		else if (sign==-1 and c!=0)
			{
			neighLocation = locNum-1;
			}
		}
	else if (c==0 and sign==-1)
		{
		neighLocation = -1;
		}
	else if (c==(xNx-1) and sign==1)
		{
		neighLocation = -1; //this is the result if there are no neighbours for the given values of the argument
		}
	else
		{
		neighLocation = locNum+sign*(int)xNt;
		}
	return neighLocation;
	}
	
long int (*neigh) (const lint& locNum, const unsigned int& direction, const signed int& sign, const unsigned int& xNt, const unsigned int & xNx);

//dt type functions
comp dtFn (const unsigned int& time)
	{
	return (time<(NT-1)? simpleTime(time+1)-simpleTime(time): 0.0);
	}
	
comp DtFn (const unsigned int& time)
	{
	if (time==(NT-1)) return (simpleTime(NT-1)-simpleTime(NT-2))/2.0;
	else if (time==0) return (simpleTime(1)-simpleTime(0))/2.0;
	else return (simpleTime(time+1)-simpleTime(time-1))/2.0;
	}
	
double dxFn (const unsigned int& space)
	{
	if (pot[0]!='3') return a;
	return (space<(N-1)? a: 0.0);
	}
	
double DxFn (const unsigned int& space)
	{
	if (pot[0]!='3') return a;
	return ((space==(N-1) || space==0)? a/2.0: a);
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//interpolate vector function
vec interpolate(vec vec_old, const unsigned int & Nt_old, const unsigned int & N_old, const unsigned int & Nt_new, const unsigned int & N_new)
	{
	unsigned int old_size = vec_old.size();
	if (old_size<2*N_old*Nt_old) {cout << "interpolate error, vec_old.size() = " << old_size << " , 2*N_old*Nt_old = " << 2*N_old*Nt_old << endl;}
	unsigned int zero_modes = old_size - 2*N_old*Nt_old;
	vec vec_new (2*N_new*Nt_new+zero_modes);
	
	unsigned int x_new, t_new, x_old, t_old;
	double exact_x_old, exact_t_old, rem_x_old, rem_t_old;
	unsigned int pos;
	int neigh_t, neigh_x, neigh_tx;	
	
	for (unsigned int l=0;l<N_new*Nt_new;l++)
		{
		t_new = intCoord(l,0,Nt_new);
		x_new = intCoord(l,1,Nt_new);
		exact_t_old = t_new*(Nt_old-1.0)/(Nt_new-1.0);
		exact_x_old = x_new*(N_old-1.0)/(N_new-1.0);
		t_old = (unsigned int)exact_t_old;
		x_old = (unsigned int)exact_x_old;
		rem_t_old = exact_t_old;
		rem_t_old -= (double)(t_old);
		rem_x_old = exact_x_old;
		rem_x_old -= (double)(x_old);
		pos = t_old + Nt_old*x_old;
		neigh_t = neigh(pos,0,1,Nt_old,N_old);
		neigh_x = neigh(pos,1,1,Nt_old,N_old);
		neigh_tx = neigh(neigh_t,1,1,Nt_old,N_old);
		if  (t_old<(Nt_old-1) )
			{
			if (neigh_t!=-1 && neigh_x)
				{
				vec_new(2*l) = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old(2*pos) \
							+ (1.0-rem_t_old)*rem_x_old*vec_old(2*neigh_x) \
							+ rem_t_old*(1.0-rem_x_old)*vec_old(2*neigh_t) \
							+ rem_t_old*rem_x_old*vec_old(2*neigh_tx);
				vec_new(2*l+1) = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old(2*pos+1)\
			 				+ (1.0-rem_t_old)*rem_x_old*vec_old(2*neigh_x+1) \
							+ rem_t_old*(1.0-rem_x_old)*vec_old(2*neigh_t+1)\
						 	+ rem_t_old*rem_x_old*vec_old(2*neigh_tx+1);
				}
			else if (neigh_t!=-1)
				{
				vec_new(2*l) = (1.0-rem_t_old)*vec_old(2*pos) \
							+ rem_t_old*vec_old(2*neigh_t);
				vec_new(2*l+1) = (1.0-rem_t_old)*vec_old(2*pos+1)\
							+ rem_t_old*vec_old(2*neigh_t+1);
				}
			else if(neigh_x!=-1)
				{
				vec_new(2*l) = (1.0-rem_x_old)*vec_old(2*pos) \
							+ rem_x_old*vec_old(2*neigh_x);
				vec_new(2*l+1) = (1.0-rem_x_old)*vec_old(2*pos+1)\
							+ rem_x_old*vec_old(2*neigh_x+1);
				}
			else
				{
				vec_new(2*l) = vec_old(2*pos);
				vec_new(2*l+1) = vec_old(2*pos+1);				
				}
			}
		}
	for (unsigned int l=0; l<zero_modes; l++)
		{
		vec_new(2*N_new*Nt_new+l) = vec_old(2*N_old*Nt_old+l);
		}
	return vec_new;
	}
	
//interpolate real vector function
vec interpolateReal(vec vec_old, const unsigned int & Nt_old, const unsigned int & N_old, const unsigned int & Nt_new, const unsigned int & N_new)
	{
	unsigned int old_size = vec_old.size();
	if (old_size<N_old*Nt_old) {cout << "interpolate error, vec_old.size() = " << old_size << " , N_old*Nt_old = " << N_old*Nt_old << endl;}
	vec vec_new(N_old*Nt_old);
	unsigned int x_new, t_new, x_old, t_old;
	double exact_x_old, exact_t_old, rem_x_old, rem_t_old;
	unsigned int pos;
	int neigh_t, neigh_x, neigh_tx;	
	
	for (unsigned int l=0;l<N_new*Nt_new;l++)
		{
		t_new = intCoord(l,0,Nt_new);
		x_new = intCoord(l,1,Nt_new);
		exact_t_old = t_new*(Nt_old-1.0)/(Nt_new-1.0);
		exact_x_old = x_new*(N_old-1.0)/(N_new-1.0);
		t_old = (unsigned int)exact_t_old;
		x_old = (unsigned int)exact_x_old;
		rem_t_old = exact_t_old;
		rem_t_old -= (double)(t_old);
		rem_x_old = exact_x_old;
		rem_x_old -= (double)(x_old);
		pos = t_old + Nt_old*x_old;
		neigh_t = neigh(pos,0,1,Nt_old,N_old);
		neigh_x = neigh(pos,1,1,Nt_old,N_old);
		neigh_tx = neigh(neigh_t,1,1,Nt_old,N_old);
		if  (t_old<(Nt_old-1) )
			{
			if (neigh_x!=-1 && neigh_t!=-1)
				{
				vec_new(l) = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old(pos) \
							+ (1.0-rem_t_old)*rem_x_old*vec_old(neigh_x) \
							+ rem_t_old*(1.0-rem_x_old)*vec_old(neigh_t) \
							+ rem_t_old*rem_x_old*vec_old(neigh_tx);
				}
			else if (neigh_t!=-1)
				{
				vec_new(l) = (1.0-rem_t_old)*vec_old(pos) \
							+ rem_t_old*vec_old(neigh_t);
				}
			else if (neigh_x!=-1)
				{
				vec_new(l) = (1.0-rem_x_old)*vec_old(pos) \
							+ rem_x_old*vec_old(neigh_x);
				}
			else
				{
				vec_new(l) = vec_old(pos);
				}
			}
		}
	return vec_new;
	}
	
//same function but 1d	
vec interpolate1d(vec vec_old, const unsigned int & N_old, const unsigned int & N_new)
	{
	unsigned int old_size = vec_old.size();
	if (old_size<N_old) cerr << "interpolate error, vec_old.size() = " << old_size << " , N_old = " << N_old << endl;
	vec vec_new (N_new);
	
	unsigned int x_old;
	double exact_x_old, rem_x_old;
	
	for (unsigned int l=0;l<N_new;l++)
		{
		exact_x_old = l*(N_old-1.0)/(N_new-1.0);
		x_old = (unsigned int)exact_x_old;
		rem_x_old = exact_x_old - (double)(x_old);
		if  (x_old<(N_old-1)) vec_new[l] = (1.0-rem_x_old)*vec_old[x_old] + rem_x_old*vec_old[x_old+1];
		else vec_new[l] = vec_old[x_old];
		}
	return vec_new;
	}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//printing functions

//print main parameters to terminal
void printParameters()
	{
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","inP","N","Na","Nb","Nc","L","Ta","Tb","Tc","R","dE","theta","reg", "epsilon");
	printf("%8s%8i%8i%8i%8i%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g\n",inP.c_str(),N,Na,Nb,Nc,L,Ta,Tb,Tc,R,dE,theta,reg,epsilon);
	printf("\n");
	}	

void printMoreParameters()
	{
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","inP","N","Na","Nb","Nc","L","Tb","R","dE","epsilon","theta","reg");
	printf("%8s%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g%8g\n",inP.c_str(),N,Na,Nb,Nc,L,Tb,R,dE,epsilon,theta,reg);
	printf("%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n","NT","Gamma","a","b","Ta","Tc","minima[0]","minima[1]","S1","mass2");
	printf("%12i%12g%12g%12g%12g%12g%12g%12g%12g%12g\n",NT,Gamma,a,b,Ta,Tc,minima[0],minima[1],S1,mass2);
	printf("\n");
	}

//simply print a real vector
void simplePrintVector(const string& printFile, vec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F.precision(16);
	F << left;
	unsigned int length = vecToPrint.size();
	for (unsigned int j=0; j<length; j++)
		{
		F << setw(25) << vecToPrint(j) << endl;
		}
	F.close();
	}

//simply print a complex vector
void simplePrintCVector(const string& printFile, cVec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F.precision(16);
	F << left;
	unsigned int length = vecToPrint.size();
	for (unsigned int j=0; j<length; j++)
		{
		F << setw(25) << real(vecToPrint(j)) << setw(25) << imag(vecToPrint(j)) << endl;
		}
	F.close();
	}
	
//print vector from time path B to file
void printVectorB (const string& printFile, vec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	unsigned int x0 = intCoord(0,1,Nb);
	F.precision(16);
	for (unsigned long int j=0; j<N*Nb; j++)
		{
		unsigned int x = intCoord(j,1,Nb);
		if (x!=x0) //this is put in for gnuplot
			{
			F << endl;
			x0 = x;
			}
		F << left;
		F << setw(24) << real(coordB(j,0)) << setw(25) << imag(coordB(j,0));
		F << setw(25) << real(coordB(j,1));
		if (vecToPrint.size()>N*Nb)
			{
			F << setw(25) << vecToPrint(2*j) << setw(25) << vecToPrint(2*j+1)  << endl;
			}
		else
			{
			F << setw(25) << vecToPrint(j) << endl;
			}
		}
	if (vecToPrint.size()>2*N*Nb)
		{
		F << endl;
		for (unsigned int k=0; k<(vecToPrint.size()-2*N*Nb);k++)
			{
			F << setw(25) << vecToPrint(2*N*Nb+k) << endl;
			}
		}
	F.close();
	}
	
//print vector to file
void printVector (const string& printFile, vec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	unsigned int x0 = intCoord(0,1,NT);
	F.precision(16);
	for (unsigned long int j=0; j<N*NT; j++)
		{
		unsigned int x = intCoord(j,1,NT);
		if (x!=x0) //this is put in for gnuplot
			{
			F << endl;
			x0 = x;
			}
		F << left;
		F << setw(25) << real(coord(j,0)) << setw(25) << imag(coord(j,0)); //note using coord for full time contour
		F << setw(25) << real(coord(j,1));
		if (vecToPrint.size()>N*NT)
			{
			F << setw(25) << vecToPrint(2*j) << setw(25) << vecToPrint(2*j+1)  << endl;
			}
		else
			{
			F << setw(25) << vecToPrint(j) << endl;
			}
		}
	if (vecToPrint.size()>2*N*NT)
		{
		F << endl;
		for (unsigned int j=0; j<(vecToPrint.size()-2*N*NT);j++)
			{
			F << setw(25) << vecToPrint(2*N*NT+j) << endl;
			}
		}
	F.close();
	}
	
//print vector to file with maximum size 300 by 300
void printReducedVector (const string& printFile, vec vecToPrint)
	{
	uint Nx = 100, Nt = 100;
	vec reducedVec = interpolate(vecToPrint,NT,N,Nt,Nx);
	ofstream F;
	F.open((printFile).c_str());
	uint x0 = intCoord(0,1,Nt);
	F.precision(16);
	for (unsigned long int j=0; j<Nx*Nt; j++)
		{
		unsigned int x = intCoord(j,1,Nt);
		if (x!=x0) //this is put in for gnuplot
			{
			F << endl;
			x0 = x;
			}
		F << left;
		F << setw(25) << real(coord(j,0)) << setw(25) << imag(coord(j,0)); //note using coord for full time contour
		F << setw(25) << real(coord(j,1));
		if (vecToPrint.size()>Nx*Nt)
			{
			F << setw(25) << vecToPrint(2*j) << setw(25) << vecToPrint(2*j+1)  << endl;
			}
		else
			{
			F << setw(25) << vecToPrint(j) << endl;
			}
		}
	if (vecToPrint.size()>2*Nx*Nt)
		{
		F << endl;
		for (unsigned int j=0; j<(vecToPrint.size()-2*Nx*Nt);j++)
			{
			F << setw(25) << vecToPrint(2*Nx*Nt+j) << endl;
			}
		}
	F.close();
	}
	
//print sparse matrix to file
void printSpmat (const string & printFile, spMat spmatToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F << left;
	F.precision(16);
	for (int l=0; l<spmatToPrint.outerSize(); ++l)
		{
		for (Eigen::SparseMatrix<double>::InnerIterator it(spmatToPrint,l); it; ++it)
			{
			F << setw(25) << it.row()+1 << setw(25) << it.col()+1 << setw(25) << it.value() << endl;
			}
		}
	F.close();
	}
	
//print p via gnuplot, using repi or pi or some such like
void gp(const string & readFile, const string & gnuplotFile) 
	{
	string prefix = "gnuplot -e \"f='";
	string middle = "'\" ";
	string suffix = " -persistent";
	string commandStr = prefix + readFile + middle + gnuplotFile + suffix;
	const char * command = commandStr.c_str();
	FILE * gnuplotPipe = popen (command,"w");
	fprintf(gnuplotPipe, "%s \n", " ");
	pclose(gnuplotPipe);
	}
	
//print simple 1:2 using gnuplot
void gpSimple(const string & readFile) 
	{
	string commandOpenStr = "gnuplot -persistent";
	const char * commandOpen = commandOpenStr.c_str();
	FILE * gnuplotPipe = popen (commandOpen,"w");
	string command1Str = "plot \"" + readFile + "\" using 1 with linespoints";
	string command2Str = "pause -1";
	const char * command1 = command1Str.c_str();
	const char * command2 = command2Str.c_str();
	fprintf(gnuplotPipe, "%s \n", command1);
	fprintf(gnuplotPipe, "%s \n", command2);
	pclose(gnuplotPipe);
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//loading functions

//load simple vector from file
vec loadSimpleVector (const string& loadFile)
	{
	unsigned int fileLength = countLines(loadFile);
	vec outputVec(fileLength);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line;
	unsigned int j=0;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			istringstream ss(line);
			ss >> outputVec(j);
			j++;
			}
		}
	F.close();
	return outputVec;
	}

//load vector from file
vec loadVector (const string& loadFile, const unsigned int& Nt, const unsigned int& Nx, const unsigned int zeroModes)
	{
	vec outputVec(2*Nt*Nx+zeroModes);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line;
	unsigned int j = 0;
	unsigned int k = 0;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			if (j<Nt*Nx)
				{
				double temp;
				istringstream ss(line);
				ss >> temp >> temp >> temp;
				ss >> outputVec(2*j) >> outputVec(2*j+1);
				j++;
				}
			else
				{
				//istringstream ss(line);
				//ss >> outputVec(2*j+k);
				outputVec(2*j+k) = 0.5; //ignoring information from previous zero mode
				k++;
				}
			}
		}
	if (j==Nt*Nx && k<zeroModes)
		{
		while(k<zeroModes)
			{
			outputVec(2*j+k) = 0.5; //obviously a random guess
			k++;
			}
		}
	if ((j+k)!=(Nt*Nx+zeroModes))
		{
		cout << "loadVector error in: " << loadFile << endl;
		cout << "j+k = " << j+k << endl;
		}
	F.close();
	return outputVec;
	}
	
//load simple vector from file
vec loadVectorColumn (const string& loadFile, const unsigned int column)
	{
	unsigned int fileLength = countLines(loadFile);
	vec outputVec(fileLength);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line, temp;
	unsigned int j=0;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			istringstream ss(line);
			for (unsigned int l=0; l<column; l++) ss >> temp;
			ss >> outputVec[j];
			j++;
			}
		}
	F.close();
	return outputVec;
	}
	
//load DDS from file
spMat loadSpmat (const string & loadFile, Eigen::VectorXi to_reserve)
	{
	unsigned int length = to_reserve.size();
	spMat M(length,length);
	M.setZero(); //just making sure
	M.reserve(to_reserve);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line;
	unsigned int nnz = 0;
	unsigned int row;
	unsigned int column;
	double value;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			istringstream ss(line);
			ss >> row >> column >> value;
			if (absolute(value)>1.0e-16)
				{
				M.insert(row-1,column-1) = value;
				nnz++;
				}
			}
		}
	if (nnz==0)
		{
		cout << "loadSpMat failed, no data in file: " << loadFile << endl;
		}
	M.makeCompressed();
	return M;
	}

//load DDS from file
unsigned int countLines(const string & file_to_count); //here so it can be used, defined in files.h
spMat cleverLoadSpmat (const string & loadFile)
	{
	unsigned int fileLength = countLines(loadFile);
	Eigen::VectorXi to_reserve(fileLength); //an overestimate
	to_reserve.setZero(fileLength);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line;
	unsigned int nnz = 0, length = 0, count = 0, row, column;
	double value;	
	Eigen::VectorXi rowVec(fileLength), columnVec(fileLength);
	vec valueVec(fileLength);
	while (getline(F, line))
		{
		if (!line.empty())
			{
			istringstream ss(line);
			ss >> row >> column >> value;
			if (absolute(value)>1.0e-16)
				{
				rowVec(count) = row-1;
				columnVec(count) = column-1;
				valueVec(count) = value;
				to_reserve(row-1) = to_reserve(row-1) + 1;
				nnz++;
				count++;
				if (row>length)
					{
					length = row;
					}
				}
			}
		}
	if (nnz==0)
		{
		cout << "loadSpMat failed, no data in file: " << loadFile << endl;
		}
	to_reserve.conservativeResize(length);
	rowVec.conservativeResize(count);
	columnVec.conservativeResize(count);
	valueVec.conservativeResize(count);
	spMat M(length,length);
	M.setZero(); //just making sure
	M.reserve(to_reserve);
	for (unsigned int l=0;l<count;l++)
		{
		M.insert(rowVec(l),columnVec(l)) = valueVec(l);
		}
	M.makeCompressed();
	return M;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//askQuestions and changeParameters

//asks initial questions to get user inputs - now defunct
void askQuestions (aqStruct & aqx )
	{
	cout << "number of loops: ";
	cin >> aqx.totalLoops;
	cout << endl;
	if (aqx.totalLoops !=1)
		{
		cout << "loop (main)parameter (N,Na,Nb,Nc,L,Tb,R,mass,lambda) ";
		cin >> aqx.loopChoice;
		cout << endl;
		cout << "min value: ";
		cin >> aqx.minValue;
		cout << endl;
		cout << "max value: ";
		cin >> aqx.maxValue;
		cout << endl;
		}
	cout << "print early (m,v,p,a,e,n)? ";
	cin >> aqx.printChoice;
	cout << endl;
	string temp = aqx.printChoice;
	if (temp.compare("n")!=0)
		{
		cout << "which run to print (0 for every run)? ";
		cin >> aqx.printRun;
		}
	}

//changes parameters according to user inputs (integer)
void changeInt (const string & parameterLabel, const int & newParameter)
	{
	if ( parameterLabel.compare("N")==0)
		{
		//Na = (int)newParameter*Na/N;
		//Nb = (int)newParameter*Nb/N;
		//Nc = (int)newParameter*Nc/N;
		//NT = Na + Nb + Nc;
		N = newParameter;
		a = L/(N-1);
		//b = Tb/(Nb-1);
		//Ta = b*(Na-1.0);
		//Tc = b*(Nc-1.0);
		}
	else if ( parameterLabel.compare("Na")==0)
		{
		Na = newParameter;
		NT = Na + Nb + Nc;
		Ta = b*(Na-1.0);
		}
	else if ( parameterLabel.compare("Nb")==0)
		{
		Nb = newParameter;
		NT = Na + Nb + Nc;
		b = Tb/(Nb-1.0);
		Ta = b*(Na-1.0);
		Tc = b*(Nc-1.0);
		}
	else if ( parameterLabel.compare("Nc")==0)
		{
		Nc = newParameter;
		NT = Nc + Nb + Nc;
		Tc = b*(Nc-1.0);
		}
	else
		{
		cout << "changeInt error" << endl;
		}
	}

//changes parameters according to user inputs (double)
void changeDouble (const string & parameterLabel, const double & newParameter)
	{
	if ( parameterLabel.compare("L")==0) //this does not changes the physics but simply the size of the box in space
		{
		L = newParameter;
		a = L/(N-1);
		}
	else if ( parameterLabel.compare("Tb")==0) //this paramter changes the physics for the periodic instanton,											//as Tb/R changes where R = R(epsilon)
		{
		b = b*newParameter/Tb;
		Ta = Ta*newParameter/Tb;
		Tc = Tc*newParameter/Tb;
		Tb = newParameter;
		if (Tb<R && pot[0]!='3'){
			angle = asin(Tb/R);
			if (2.0*(1.5*Tb*tan(angle))<L) L=2.0*(1.5*Tb*tan(angle));
			a = L/(N-1.0);
			}
		}
	else if ( parameterLabel.compare("R")==0) //this parameter changes the initial guess
		{
		L = L*newParameter/R; //all length scales scale with R
		a = a*newParameter/R;
		b = b*newParameter/R;
		Ta = Ta*newParameter/R;
		Tb = Tb*newParameter/R;
		Tc = Tc*newParameter/R;
		R = newParameter;
		}
	else if ( parameterLabel.compare("dE")==0) //this parameter changes the physics of the potential
													//but it does not change Tb/R, where R(epsilon)
		{
		R = R*dE/newParameter; //R scales with 1/dE and the other length scales scale with R
		L = L*dE/newParameter;
		a = a*dE/newParameter;
		b = b*dE/newParameter;
		Ta = Ta*dE/newParameter;
		Tb = Tb*dE/newParameter;
		Tc = Tc*dE/newParameter;
		epsilon = epsilon*newParameter/dE;
		dE = newParameter;
		}
	else
		{
		cout << "changeDouble error" << endl;
		}
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//vector manipulation functions

//complexify a real vector
cVec vecComplex(vec realVec, const unsigned int & tDim)
	{
	cVec complexVec(tDim);
	if (realVec.size() >= (2*tDim) && realVec.size() < (2*tDim+3))
		{
		for (unsigned int l=0; l<tDim; l++)
			{
			complexVec(l) = realVec(2*l) + ii*realVec(2*l+1);
			}
		}
	else
		{
		cout << "vecComplex error";
		}
	return complexVec;
	}
	
//make a complex vector real
vec vecReal(cVec complexVec, const unsigned int &  tDim)
	{
	vec realVec(2*tDim);
	if (complexVec.size() == tDim)
		{
		for (unsigned int l=0; l<tDim; l++)
			{
			realVec(2*l) = real(complexVec(l));
			realVec(2*l+1) = imag(complexVec(l));
			}
		}
	else
		{
		cout << "vecReal error";
		}
	return realVec;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//fourier transform type functions

//h the matrix from dl[7]
mat hFn(const unsigned int & xN, const double & xa, const double & xmass2, const string& potential=pot)
	{
	mat xh(xN,xN);	xh = Eigen::MatrixXd::Zero(xN,xN);
	double diag = xmass2 + 2.0/pow(xa,2.0);
	double offDiag1 = -1.0/pow(xa,2.0);
	double offDiag2;
	if (pot[0]=='3') offDiag2 = -pow(2.0,0.5)/pow(xa,2.0);	
	else			 offDiag2 = -1.0/pow(xa,2.0);
	for (unsigned int l=0; l<xN; l++)
		{
		if (l==0)
			{
			xh(l,l) = diag;
			xh(l,l+1) = offDiag2;			
			}
		else if (l==(xN-1))
			{
			xh(l,l) = diag;
			xh(l,l-1) = offDiag2;
			}
		else
			{
			xh(l,l) = diag;
			if ((l+1)==(xN-1))	xh(l,l+1) = offDiag2;
			else				xh(l,l+1) = offDiag1;
			if ((l-1)==0)		xh(l,l-1) = offDiag2;
			else				xh(l,l-1) = offDiag1;
			}
		}
	return xh;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//misc functions


