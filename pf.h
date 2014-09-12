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
unsigned int N;
unsigned int Na;
unsigned int Nb;
unsigned int Nc;
double dE;
double theta;

//derived parameters
unsigned int NT;
double epsilon;
double R; //size of bubble
double Tb;
double angle; //not a primary parameter, just used to make L
double L;
double a; //step sizes in each spatial dimension
double b; //step sizes in time
double Ta;
double Tc;

//determining number of runs
double closenessA; //action
double closenessS; //solution (i.e. minusDS)
double closenessD; //delta
double closenessC; //calculation

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

//potential related functions
comp V (const comp & phi)
	{
	comp Vphi = pow(pow(phi,2)-1.0,2)/8.0 - epsilon*(phi-1.0)/2.0;
	return Vphi;
	}
	
//potential with degenerate minima
comp Va (const comp & phi)
	{
	comp Vaphi = pow(pow(phi,2)-1.0,2)/8.0;
	return Vaphi;
	}
	
//change to degenerate potential
comp Vb (const comp & phi)
	{
	comp Vbphi = -epsilon*(phi-1.0)/2.0;
	return Vbphi;
	}

//first derivative of V
comp dV (const comp & phi)
	{
	comp dVphi = phi*(pow(phi,2)-1.0)/2.0 - epsilon/2.0;
	return dVphi;
	}
	
//second derivative of V
comp ddV (const comp & phi)
	{
	comp ddVphi = (3.0*pow(phi,2)-1.0)/2.0;
	return ddVphi;
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
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//printing functions

//print main parameters to terminal
void printParameters()
	{
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","N","Na","Nb","Nc","L","Tb","R","dE","epsilon","theta");
	printf("%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g\n",N,Na,Nb,Nc,L,Tb,R,dE,epsilon,theta);
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
			F << setw(25) << real(coordB(j,r)) << setw(25) << imag(coordB(j,r));
			}
		if (vecToPrint.size()>N*Nb)
			{
			F << setw(25) << vecToPrint(2*j) << setw(25) << vecToPrint(2*j+1)  << endl;
			}
		else
			{
			F << setw(25) << vecToPrint(j) << endl;
			}
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
	string inputChoice;
	unsigned int fileNo;
	double maxTheta;
	unsigned int totalLoops;
	string loopChoice;
	double minValue;
	double maxValue;
	string printChoice;
	unsigned int printRun;
	};

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
	else if ( parameterLabel.compare("Tb")==0) //this paramter changes the physics for the periodic instanton,
												//as Tb/R changes where R = R(epsilon)
		{
		b = b*newParameter/Tb;
		Ta = Ta*newParameter/Tb;
		Tc = Tc*newParameter/Tb;
		Tb = newParameter;
		angle = asin(Tb/R);
		L = 3*R;
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
	else if ( parameterLabel.compare("epsilon")==0) //this parameter changes the physics of the potential
													//but it does not change Tb/R, where R(epsilon)
		{
		R = R*epsilon/newParameter; //R scales with 1/epsilon and the other length scales scale with R
		L = L*epsilon/newParameter;
		a = a*epsilon/newParameter;
		b = b*epsilon/newParameter;
		Ta = Ta*epsilon/newParameter;
		Tb = Tb*epsilon/newParameter;
		Tc = Tc*epsilon/newParameter;
		epsilon = newParameter;
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
