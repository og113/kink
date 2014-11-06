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
#include "gnuplot_i.hpp"

using namespace std;

typedef unsigned long int lint;
typedef complex<double> comp;
typedef vector<unsigned int> intVec;
typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;
typedef Eigen::Triplet<double> triplet;

complex<double> i(0.0,1.0);
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
vector<double> root;
double mass2; //as derived from V''

//determining number of runs
double closenessA; //action
double closenessS; //solution (i.e. minusDS)
double closenessSM; //solution max
double closenessD; //delta
double closenessC; //calculation
double closenessE; //energy change
double closenessL; //linearisation of energy
double closenessT; //true energy versus linear energy
double closenessP; //checking lattice small enough for momenta
double closenessR; //regularization

//parameters determining input phi
//struct to hold answers to questions
struct aqStruct
	{
	string inputChoice;
	string inputFile;
	double maxTheta;
	unsigned int totalLoops;
	string loopChoice;
	double minValue;
	double maxValue;
	string printChoice;
	unsigned int printRun;
	};
aqStruct aq; //struct to hold user responses
string inP; //b for bubble, p for periodic instaton, f for from file
string pot; //pot[0] gives 1 or 2, pot[1] gives r (regularised) or n (not)
double A; //gives parameter in V2, equal to 0.4 in DL
double reg; //small parameter multiplying regulatory term
string inF; //file input from where, m for main, p for pi
int firstLoop, lastLoop; //first and last loop to load from, not unsigned for comparison later
double alpha; //gives span over which tanh is used
double open; //value of 0 assigns all weight to boundary, value of 1 to neighbour of boundary
double amp; //ammount of negative eigenvector added to bubble for Tb>R
double negVal; //the negative eigenvalue
unsigned int negEigDone; //has the negEig been found before? 1 if yes, 0 if no
string zmt; //dictates how time zero mode is dealt with
string zmx; //dictates how x zero mode is dealt with
double epsilon0; //the value of epsilon when the minima are degenerate

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//simple generic functions

//factorial function
int factorial(const int& f_input)
	{
	int f_result = 1;
	for (int l=0; l<f_input; l++)
		{
		f_result *= (l+1);
		}
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
	
//function giving location of smallest element of a unsigned int vector
unsigned int smallestFn(const vector <unsigned int> & inVector)
	{
	unsigned int loc = 0;
	for(unsigned int l=1;l<inVector.size();l++)
		{
		if (inVector[l]<inVector[loc])
			{
			loc = l;
			}
		}
	return loc;
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
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//potential related functions
comp V1_params (const comp & phi, const double & epsi, const double & aa) //aa is intentionally unused
	{
	return pow(pow(phi,2)-1.0,2.0)/8.0 - epsi*(phi-1.0)/2.0;
	}

comp V1 (const comp & phi)
	{
	return V1_params(phi,epsilon,0.0);
	}
	
comp Z (const comp & phi)
	{
	return exp(-pow(phi,2.0))*(phi + pow(phi,3.0) + pow(phi,5.0));
	}
	
comp V2_params (const comp & phi, const double & epsi, const double & aa)
	{
	return 0.5*pow(phi+1.0,2.0)*(1.0-epsi*Z((phi-1.0)/aa));
	}
	
comp V2 (const comp & phi)
	{
	return V2_params(phi,epsilon,A);
	}
	
comp VrFn (const comp & phi, const double & minimaL, const double & minimaR)
	{
	return pow(phi-minimaL,4.0)*pow(phi-minimaR,4.0)/4.0;
	}
	
comp (*V) (const comp & phi); //a function pointer

comp (*V_params) (const comp & phi, const double & epsi, const double & aa);

//potentials with degenerate minima, and without regularisation
comp V10 (const comp & phi)
	{
	return V1_params(phi,epsilon0,0.0);
	}
	
comp V20 (const comp & phi)
	{
	return V2_params(phi,epsilon0,A); //epsilon0 needs checking
	}
	
comp (*V0) (const comp & phi);
	
//change to degenerate potential, still without regulatisation
comp V1e (const comp & phi)
	{
	return -epsilon*(phi-1.0)/2.0;
	}

comp V2e (const comp & phi)
	{
	return V2(phi)-V20(phi);
	}
	
comp (*Ve) (const comp & phi);
	
	
//first derivative of Vs
comp dV1_params(const comp & phi, const double & epsi, const double & aa) //aa is intentionally unused
	{
	return phi*(pow(phi,2)-1.0)/2.0 - epsi/2.0;
	}

comp dV1 (const comp & phi)
	{
	return dV1_params(phi,epsilon,0.0);
	}

comp dZ (const comp & phi)
	{
	return exp(-pow(phi,2.0))*(1.0 + pow(phi,2.0) + 3.0*pow(phi,4.0) -2.0*pow(phi,6.0));
	}
	
comp dV2_params(const comp & phi, const double & epsi, const double & aa)
	{
	return (phi+1.0)*(1.0-epsi*Z((phi-1.0)/aa)) - 0.5*pow(phi+1.0,2.0)*(epsi/aa)*dZ((phi-1.0)/aa);
	}

comp dV2 (const comp & phi)
	{
	return dV2_params(phi,epsilon,A);
	}
	
comp dVrFn (const comp & phi, const double & minimaL, const double & minimaR)
	{
	return pow(phi-minimaL,3.0)*pow(phi-minimaR,4.0) + pow(phi-minimaL,4.0)*pow(phi-minimaR,3.0);
	}
	
comp (*dV) (const comp & phi);

comp (*dV_params)(const comp & phi, const double & epsi, const double & aa);
	
//second derivative of V
comp ddV1_params (const comp & phi, const double & epsi, const double & aa)
	{ 
	return (3.0*pow(phi,2)-1.0)/2.0;
	}

comp ddV1 (const comp & phi)
	{ 
	return ddV1_params(phi,epsilon,A);
	}

comp ddZ (const comp & phi)
	{
	return exp(-pow(phi,2.0))*2.0*pow(phi,3.0)*(5.0 - 9.0*pow(phi,2.0) + 2.0*pow(phi,4.0));
	}

comp ddV2_params (const comp & phi, const double & epsi, const double & aa)
	{ 
	return (1.0-epsi*Z((phi-1.0)/aa)) - (phi+1.0)*(epsi/aa)*dZ((phi-1.0)/aa)\
					+ 0.5*pow(phi+1.0,2.0)*(epsi/pow(aa,2.0))*ddZ((phi-1.0)/aa);
	}
	
comp ddV2 (const comp & phi)
	{ 
	return ddV2_params(phi,epsilon,A);
	}	
	
comp ddVrFn (const comp & phi, const double & minimaL, const double & minimaR)
	{
	return 3.0*pow(phi-minimaL,2.0)*pow(phi-minimaR,4.0) + 8.0*pow(phi-minimaL,3.0)*pow(phi-minimaR,3.0)\
				+ 3.0*pow(phi-minimaL,4.0)*pow(phi-minimaR,2.0);
	}

comp (*ddV) (const comp & phi);

comp (*ddV_params)(const comp & phi, const double & epsi, const double & aa);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//functions for calculating roots
	
//dV as gsl ready functions
struct f_gsl_params { double epsi; double aa;};

double f_gsl (double x, void * parameters) 
	{
	struct f_gsl_params * params = (struct f_gsl_params *)parameters;
	double epsi = (params->epsi);
	double aa = (params->aa);
	return real(dV_params(x,epsi,aa));
	}
	
double df_gsl (double x, void * parameters) 
	{
	struct f_gsl_params * params = (struct f_gsl_params *)parameters;
	double epsi = (params->epsi);
	double aa = (params->aa);
	return real(ddV_params(x,epsi,aa));
	}
	
void fdf_gsl (double x, void * parameters, double * f, double* df) 
	{
	struct f_gsl_params * params = (struct f_gsl_params *)parameters;
	double epsi = (params->epsi);
	double aa = (params->aa);
	*f =  real(dV_params(x,epsi,aa));
	*df = real(ddV_params(x,epsi,aa));
	}
	
//energy change gsl functions : V(root[2])-V(root[0])-dE
struct ec_gsl_params {double aa; double root0; double root2; double de; };

double ec_gsl (double epsi, void * parameters)
	{
	struct ec_gsl_params * params = (struct ec_gsl_params *)parameters;
	double aa = (params->aa);
	double root2 = (params->root2);
	double root0 = (params->root0);
	double de = (params->de);
	return real(V_params(root0,epsi,aa) - V_params(root2,epsi,aa) - de);
	}
	
//dV0
struct void_gsl_params{};
//dV0 as a gsl function
double dV0_gsl (double x, void * parameters)
	{
	return real(dV_params(x,epsilon0,A));
	}
	
//ddV0 as a gsl_function
double ddV0_gsl (double x, void * parameters)
	{
	return real(ddV_params(x,epsilon0,A));
	}
	
//dV0ddV0 as a void gsl_fdf
void dV0ddV0_gsl (double x, void * parameters, double * f, double* df) 
	{
	*f =  real(dV_params(x,epsilon0,A));
	*df = real(ddV_params(x,epsilon0,A));
	}
	
//S1 integrand
double s1_gsl (double x, void * parameters)
	{
	return pow(2.0*real(V0(x)),0.5);
	}

//rho integrand
double rho_gsl (double x, void * parameters)
	{
	return pow(2.0*real(V0(x)),-0.5);
	}
	
//function to give the three roots of FDF given lower and upper limits on them and a number of loops to try
vector <double> minimaFn (gsl_function_fdf * xFDF, const double & lowLimit, const double & highLimit,\
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
			cout << "minimaFn error: only found " << root.size() << " roots, not 3" << endl;
			}
	return roots;
	}
	
//program to find epsilon given a gsl function fdf and dE
void epsilonFn (gsl_function_fdf * xFDF, gsl_function * xEC, double * xdE, double * xEpsilon, vector<double>* xRoot)
	{
	double closenessdE =  DBL_MIN;
	vector<double> dE_test(1);	dE_test[0] = 1.0;
	double newdE = dE;
	struct f_gsl_params * Fparameters = (struct f_gsl_params *) (*xFDF).params;
	struct ec_gsl_params * ECparameters = (struct ec_gsl_params *) (*xEC).params;
	unsigned int counter = 0;
	unsigned int maxCounter = 100;
	while (dE_test.back()>closenessdE)
		{
		//find roots of ec(epsilon)=0
		*xEpsilon = brentRootFinder(xEC,*xEpsilon,*xEpsilon/2.0,*xEpsilon*2.0);
		//assign new value of epsilon to xFDF
		(*Fparameters).epsi = *xEpsilon;
		(*xFDF).params = Fparameters;
		//finding new roots of dV(phi)=0
		*xRoot = minimaFn(xFDF, -3.0, 3.0, 20);
		sort((*xRoot).begin(),(*xRoot).end());
		//assign new roots to xECDF
		(*ECparameters).root0 = (*xRoot)[0];
		(*ECparameters).root2 = (*xRoot)[2];
		(*xEC).params = ECparameters;
		//evaluating new dE
		newdE = (*(*xEC).function)(*xEpsilon,ECparameters) + dE;
		//evaluating test
		dE_test.push_back(absolute((newdE-(*xdE))/(*xdE)));
		counter++;
		//test if too many runs
		if (counter>maxCounter)
			{
			cout << "epsilonFn error, more that " << maxCounter << " loops, consider reducing closenessdE" << endl;
			cout << "dE_test.back() = " << dE_test.back() << endl;
			}
		}
	*xdE = newdE;
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//lattice functions

//function which find integer coordinates in 2d
unsigned int intCoord(const unsigned int& locNum, const int& direction, const unsigned int& xNt)
	{
	unsigned int XintCoord;
	unsigned int x = floor(locNum/xNt);
	if (direction==1)
		{
		XintCoord = x;
		}
	else
		{
		XintCoord = locNum - x*xNt;
		}
	return XintCoord;	
	}
	
//simple time
comp simpleTime (const unsigned int& time)
	{
	comp xTime;
	if ( time < Na)
		{
		double temp = (double)time;
		temp -= (double)Na;
		xTime = b*temp + i*Tb;
		}
	else if (time < (Na+Nb))
		{
		double temp = (double)time;
		temp -= Na; //as complex doesn't support (complex double)*integer (though it does support double*integer added to a complex double) - and as int to double seems to cock up here (perhaps because the integers are unsigned)
		xTime = i*Tb - i*b*temp;
		}
	else
		{
		double temp = (double)time;
		temp -= (double)Na;
		temp -= (double)Nb;
		xTime = b*(temp+1.0); //the 1.0 is because the corner is part of the vertical contour
		}
	return xTime;
	}
	
//simple space
double simpleSpace (const unsigned int& space)
	{
	double xSpace = -L/2.0 + space*a;
	return xSpace;
	}
	
//gives values of coordinates in whole spacetime
comp coord(const unsigned int& locNum,const int& direction)
	{
	comp xCoord;
	if (direction==0)
		{
		unsigned int t = intCoord(locNum,0,NT);
		xCoord = simpleTime(t);
		}
	if (direction==1)
		{
		unsigned int x = intCoord(locNum,1,NT);
		xCoord = simpleSpace(x);
		}
	return xCoord;
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
		XcoordA = b*temp + i*Tb;
		}
	if (direction==1)
		{
		unsigned int x = intCoord(locNum,1,Na);
		double temp = (double)x;
		XcoordA = -L/2.0 + temp*a;
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
		XcoordB = i*(Tb - b*temp);
		}
	if (direction==1)
		{
		unsigned int x = intCoord(locNum,1,Nb);
		double temp = (double)x;
		XcoordB = -L/2.0 + temp*a;
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
		double temp = (double)x;
		XcoordC = -L/2.0 + temp*a;
		}
	return XcoordC;
	}
	
long int neigh(const lint& locNum, const unsigned int& direction, const signed int& sign, const unsigned int& xNt) //periodic in space but not time, degree refers to the number of neighbours, 1 is for just positive neighbours, 2 is for both
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
		neighLocation = locNum+(N-1)*xNt;
		}
	else if (c==(N-1) and sign==1)
		{
		neighLocation = locNum-(N-1)*(int)xNt;
		}
	else
		{
		neighLocation = locNum+sign*(int)xNt;
		}
	return neighLocation;
	}
	
comp dtFn (const unsigned int& time)
	{
	comp xdt;
	if (time<(NT-1))
		{
		xdt = simpleTime(time+1)-simpleTime(time);
		}
	else
		{
		xdt = 0;
		}
	return xdt;
	}
	
comp DtFn (const unsigned int& time)
	{
	comp xDt;
	if (time==(NT-1))
		{
		xDt = (simpleTime(time)-simpleTime(time-1))/2.0;
		}
	else if (time==0)
		{
		xDt = (simpleTime(time+1) - simpleTime(time))/2.0;
		}
	else
		{
		xDt = (simpleTime(time+1)-simpleTime(time-1))/2.0;
		}
	return xDt;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//printing functions

//print main parameters to terminal
void printParameters()
	{
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","inP","N","Na","Nb","Nc","L","Tb","R","dE","epsilon","theta","reg");
	printf("%8s%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g%8g\n",inP.c_str(),N,Na,Nb,Nc,L,Tb,R,dE,epsilon,theta,reg);
	printf("\n");
	}	

void printMoreParameters()
	{
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","inP","N","Na","Nb","Nc","L","Tb","R","dE","epsilon","theta","reg");
	printf("%8s%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g%8g\n",inP.c_str(),N,Na,Nb,Nc,L,Tb,R,dE,epsilon,theta,reg);
	printf("%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n","NT","Gamma","a","b","Ta","Tc","root[0]","root[1]","root[2]","mass2");
	printf("%12i%12g%12g%12g%12g%12g%12g%12g%12g%12g\n",NT,Gamma,a,b,Ta,Tc,root[0],root[1],root[2],mass2);
	printf("\n");
	}
	
//print action and its constituents to the terminal
void printAction ( const comp& Kinetic, const comp& potL, const comp& potE)
	{
	comp action = Kinetic + potL + potE;
	printf("%16s%16s%16s%16s%16s%16s%16s%16s\n","re(kinetic)","im(kinetic)","re(potL)","im(potL)","re(potE)","im(potE)","re(action)","im(action)");
	printf("%16g%16g%16g%16g%16g%16g%16g%16g\n",real(Kinetic),imag(Kinetic),real(potL),imag(potL),real(potE),imag(potE),real(action),imag(action));
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
	
//print repi via gnuplot
void gpSimple(const string & readFile) 
	{
	string commandOpenStr = "gnuplot -persistent";
	const char * commandOpen = commandOpenStr.c_str();
	FILE * gnuplotPipe = popen (commandOpen,"w");
	string command1Str = "plot \"" + readFile + "\" using 1 with lines";
	string command2Str = "pause -1";
	const char * command1 = command1Str.c_str();
	const char * command2 = command2Str.c_str();
	fprintf(gnuplotPipe, "%s \n", command1);
	fprintf(gnuplotPipe, "%s \n", command2);
	pclose(gnuplotPipe);
	}
	
//print repi via gnuplot
void gpSimple2(const string & readFile) 
	{
	string commandOpenStr = "gnuplot -persistent";
	const char * commandOpen = commandOpenStr.c_str();
	FILE * gnuplotPipe = popen (commandOpen,"w");
	string command1Str = "plot \"" + readFile + "\" using 1:2 with points";
	string command2Str = "pause -1";
	const char * command1 = command1Str.c_str();
	const char * command2 = command2Str.c_str();
	fprintf(gnuplotPipe, "%s \n", command1);
	fprintf(gnuplotPipe, "%s \n", command2);
	pclose(gnuplotPipe);
	}

	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//loading functions

//load vector from file
vec loadVector (const string& loadFile, const unsigned int& Nt, const unsigned int zeroModes)
	{
	vec outputVec(2*Nt*N+zeroModes);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line;
	unsigned int j = 0;
	unsigned int k = 0;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			if (j<Nt*N)
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
	if (j==Nt*N && k==(zeroModes-1))
		{
		outputVec(2*j+k) = 0.5; //obviously a random guess
		k++;
		}
	if ((j+k)!=(Nt*N+zeroModes))
		{
		cout << "loadVector error in: << " << loadFile << endl;
		cout << "j+k = " << j+k << endl;
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
			if (absolute(value)>DBL_MIN)
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
		Na = (int)newParameter*Na/N;
		Nb = (int)newParameter*Nb/N;
		Nc = (int)newParameter*Nc/N;
		NT = Na + Nb + Nc;
		N = newParameter;
		a = L/(N-1);
		b = Tb/(Nb-1);
		Ta = b*(Na-1.0);
		Tc = b*(Nc-1.0);
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
		angle = asin(Tb/R);
		if (2.0*(1.5*Tb*tan(angle))<L) { L=2.0*(1.5*Tb*tan(angle));}
		a = L/(N-1.0);
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
			complexVec(l) = realVec(2*l) + i*realVec(2*l+1);
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
mat hFn(const unsigned int & xN, const double & xa, const double & mass2)
	{
	mat xh(xN,xN);	xh = Eigen::MatrixXd::Zero(xN,xN);
	double diag = mass2 + 2.0/pow(xa,2.0); //diagonal terms
	double offDiag = -1.0/pow(xa,2.0); //off diagonal terms
	for (unsigned int l=0; l<xN; l++)
		{
		if (l==0)
			{
			xh(l,l) = diag;
			xh(l,l+1) = offDiag;			
			}
		else if (l==(xN-1))
			{
			xh(l,l) = diag;
			xh(l,l-1) = offDiag;
			}
		else
			{
			xh(l,l) = diag;
			xh(l,l+1) = offDiag;
			xh(l,l-1) = offDiag;	
			}
		}
	return xh;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//misc functions


