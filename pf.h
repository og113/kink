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
#include <gsl/gsl_poly.h>

using namespace std;

typedef unsigned long int lint;
typedef complex<double> comp;
typedef vector<unsigned int> intVec;
typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;

complex<double> i(0.0,1.0);
#define pi 3.14159265359

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//primary parameters
unsigned int N = 80; //number of points in each spatial dimension
unsigned int Na = (int)(1.2*N);
unsigned int Nb = (int)(1.0*N);
unsigned int Nc = 2;
double R = 50.0; //size of bubble
double mass = 1.0; 
double lambda = 0.1;
double Tb = R;
double angle = asin(Tb/R); //not a primary parameter, just used to make L
double L = 3.0*R;
double Ltemp = 2.0*(1.5*Tb*tan(angle));
double theta = 0.0;

//derived quantities
unsigned int NT = Na + Nb + Nc;
double a = L/(N-1.0); //step sizes in each spatial dimension
double b = Tb/(Nb-1.0); //step sizes in time
double Ta = b*(Na-1.0);
double Tc = b*(Nc-1.0);
double v =  mass*pow(lambda,-0.5); //vacuum phi
double X = mass*R; //the large thin-wall parameter, the ratio of the size of the bubble to the size of the wall
double epsilon = 2.0*pow(mass,3)/lambda/R/3.0; //energy difference

//determining number of runs
double closenessA = pow(10,-4.0);
double closenessV = v*pow(10,-3.0);
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
	
//gives values of coordinates in whole spacetime
comp coord(const unsigned int& locNum,const int& direction)
	{
	comp Xcoord;
	if (direction==0)
		{
		unsigned int t = intCoord(locNum,0,NT);
		if ( t < Na)
			{
			double temp = (double)t;
			temp -= (double)Na;
			Xcoord = b*temp + i*Tb;
			}
		else if (intCoord(locNum,0,NT) < (Na+Nb))
			{
			double temp = (double)t;
			temp -= Na; //as complex doesn't support (complex double)*integer (though it does support double*integer added to a complex double) - and as int to double seems to cock up here (perhaps because the integers are unsigned)
			Xcoord = i*Tb - i*b*temp;
			}
		else
			{
			double temp = (double)t;
			temp -= (double)Na;
			temp -= (double)Nb;
			Xcoord = b*temp;
			}
		}
	if (direction==1)
		{
		unsigned int x = intCoord(locNum,1,NT);
		double temp = (double)x;
		Xcoord = -L/2.0 + temp*a;
		}
	return Xcoord;
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
		if (sign==-1 and c!=0)
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
		neighLocation = locNum-(N-1)*xNt;
		}
	else
		{
		neighLocation = locNum+sign*xNt;
		}
	return neighLocation;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//printing functions

//print main parameters to terminal
void printParameters()
	{
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","N","Na","Nb","Nc","L","Tb","R","mass","lambda","theta");
	printf("%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g\n",N,Na,Nb,Nc,L,Tb,R,mass,lambda,theta);
	}
	
//print action and its constituents to the terminal
void printAction ( const comp& Kinetic, const comp& potL, const comp& potE)
	{
	comp action = Kinetic + potL + potE;
	printf("%16s%16s%16s%16s%16s%16s%16s%16s\n","re(kinetic)","im(kinetic)","re(potL)","im(potL)","re(potE)","im(potE)","re(action)","im(action)");
	printf("%16g%16g%16g%16g%16g%16g%16g%16g\n",real(Kinetic),imag(Kinetic),real(potL),imag(potL),real(potE),imag(potE),real(action),imag(action));
	}
	
//print vector from time path B to file
void printVectorB (const string& printFile, vec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	double x0 = intCoord(0,1,Nb);
	F.precision(16);
	for (unsigned long int j=0; j<N*Nb; j++)
		{
		double x = intCoord(j,1,Nb);
		if (x!=x0) //this is put in for gnuplot
			{
			F << endl;
			x0 = x;
			}
		F << left;
		for (int r=0; r<2; r++)
			{
			F << setw(25) << real(coordB(j,r)) << setw(25) << imag(coordB(j,r)); //note using coord for full time contour
			}
		F << setw(25) << vecToPrint(2*j) << setw(25) << vecToPrint(2*j+1)  << endl;
		}
	F.close();
	}
	
//print vector to file
void printVector (const string& printFile, vec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	double x0 = intCoord(0,1,NT);
	F.precision(16);
	for (unsigned long int j=0; j<N*NT; j++)
		{
		double x = intCoord(j,1,NT);
		if (x!=x0) //this is put in for gnuplot
			{
			F << endl;
			x0 = x;
			}
		F << left;
		F << setw(25) << real(coord(j,0)) << setw(25) << imag(coord(j,0)); //note using coord for full time contour
		F << setw(25) << real(coord(j,1)) << setw(25) << imag(coord(j,1));
		F << setw(25) << vecToPrint(2*j) << setw(25) << vecToPrint(2*j+1)  << endl;
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
			F << setw(25) << it.row()+1 << setw(15) << it.col()+1 << setw(25) << it.value() << endl;
			}
		}
	F.close();
	}

	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//askQuestions and changeParameters

//struct to hold answers to questions
struct aqStruct
	{
	unsigned int fileNo;
	double maxTheta;
	unsigned int totalLoops;
	string loopChoice;
	double minValue;
	double maxValue;
	string printChoice;
	unsigned int printRun;
	};

//asks initial questions to get user inputs
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
	if ( parameterLabel.compare("L")==0)
		{
		L = newParameter;
		a = L/(N-1);
		}
	else if ( parameterLabel.compare("Tb")==0)
		{
		Tb = newParameter;
		b = Tb/(Nb-1.0);
		Ta = b*(Na-1.0);
		Tc = b*(Nc-1.0);
		angle = asin(Tb/R);
		L = 3*R;
		if (2.0*(1.5*Tb*tan(angle))<L) { L=2.0*(1.5*Tb*tan(angle));}
		a = L/(N-1.0);
		}
	else if ( parameterLabel.compare("R")==0)
		{
		Tb = newParameter*Tb/R;
		R = newParameter;
		X = mass*R;
		epsilon = 2.0*pow(mass,3)/lambda/R/3.0;
		L = 3*R;
		if (2.0*(1.5*Tb*tan(angle))<L) { L=2.0*(1.5*Tb*tan(angle));}
		a = L/(N-1);
		}
	else if ( parameterLabel.compare("mass")==0)
		{
		mass = newParameter;
		v =  mass*pow(lambda,-0.5);
		X = mass*R;
		epsilon = 2.0*pow(mass,3)/lambda/R/3.0;
		}
	else if ( parameterLabel.compare("lambda")==0)
		{
		lambda = newParameter;
		v =  mass*pow(lambda,-0.5);
		epsilon = 2.0*pow(mass,3)/lambda/R/3.0;
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
	if (realVec.size() == (2*tDim+1) || realVec.size() == 2*tDim)
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
//miscelaneous other functions
	
//function to check whether the solution has converged
void convergence (const int& runsCount, const double& actionTest, const double& vectorTest, const clock_t& Clock, clock_t* Wait, char* printWait, string* printChoice, unsigned int* printRun, const comp& Kinetic, const comp& potL, const comp& potE, bool boolWait)
	{
	//stopping newton-raphson loop if action doesn't converge after 1000 loops
	if (runsCount > 1000)
		{
		cout << "over 1000 runs without convergence" << endl;
		cout << "action_test = " << actionTest << ", closennesA = " << closenessA << endl;
		cout << "vector_test = " << vectorTest << ", closennesV = " << closenessV << endl;
		}

//prints runs_count if looping is taking ages
	if((Clock-*Wait)/1.0e6>600 and not(boolWait))
		{
		cout << "number of newton-raphson loops = " << runsCount << endl;
		*Wait = Clock;
		cout << "print phi and action on the next loop? (y/n)" << endl;
		cin >> *printWait;
		if (*printWait == 'y')
			{
			comp Action = Kinetic+potL+potE;
			*printChoice = 'p';
			*printRun = runsCount+1;
			cout << left;
			cout << "kinetic = " << Kinetic << endl;
			cout << "pot_lambda = " << potL << endl;
			cout << "pot_epsilon = " << potE << endl;
			cout << "action = " << Action << endl;
			}
		}
	}
	
