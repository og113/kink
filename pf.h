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
#include <gsl/gsl_poly.h>
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
double dE;
double theta;

//derived parameters
unsigned int NT;
double epsilon;
double R; //size of bubble
double Tb;  //b section includes both corner points
double Gamma; //equals exp(-theta)
double angle; //not a primary parameter, just used to make L
double L;
double a; //step sizes in each spatial dimension
double b; //step sizes in time
double Ta;
double Tc;
vector<double> root(3);

//determining number of runs
double closenessA; //action
double closenessS; //solution (i.e. minusDS)
double closenessSM; //solution max
double closenessD; //delta
double closenessC; //calculation

//parameters determining input phi
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
aqStruct aq; //struct to hold user responses
string inP; //b for bubble, p for periodic instaton, f for from file
double alpha; //gives span over which tanh is used
double open; //value of 0 assigns all weight to boundary, value of 1 to neighbour of boundary
double amp; //ammount of negative eigenvector added to bubble for Tb>R

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
		xTime = b*temp;
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
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","inP","N","Na","Nb","Nc","L","Tb","R","dE","epsilon","theta");
	printf("%8s%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g\n",inP.c_str(),N,Na,Nb,Nc,L,Tb,R,dE,epsilon,theta);
	printf("\n");
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
		for (int r=0; r<2; r++)
			{
			F << setw(20) << real(coordB(j,r)) << setw(20) << imag(coordB(j,r));
			}
		if (vecToPrint.size()>N*Nb)
			{
			F << setw(20) << vecToPrint(2*j) << setw(20) << vecToPrint(2*j+1)  << endl;
			}
		else
			{
			F << setw(20) << vecToPrint(j) << endl;
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
		F << setw(22) << real(coord(j,0)) << setw(22) << imag(coord(j,0)); //note using coord for full time contour
		F << setw(22) << real(coord(j,1)) << setw(22) << imag(coord(j,1));
		if (vecToPrint.size()>N*NT)
			{
			F << setw(22) << vecToPrint(2*j) << setw(22) << vecToPrint(2*j+1)  << endl;
			}
		else
			{
			F << setw(22) << vecToPrint(j) << endl;
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
			F << setw(22) << it.row()+1 << setw(22) << it.col()+1 << setw(22) << it.value() << endl;
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
	FILE * gnuplotPipe = popen (command,"w");
	string command1Str = "plot \"" + readFile + "\" with lines";
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
vec loadVector (const string& loadFile, const unsigned int& Nt)
	{
	vec outputVec(2*Nt*N+1);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line;
	unsigned int j = 0;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			double temp;
			istringstream ss(line);
			ss >> temp >> temp >> temp >> temp;
			ss >> outputVec(2*j) >> outputVec(2*j+1);
			j++;
			}
		}
	if (j!=Nt*N)
		{
		cout << "loadVector error" << endl;
		}
	outputVec(2*Nt*N) = 0.5; //lagrange multiplier to remove zero mode
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
			if (absolute(value)>2.0e-16)
				{
				M.insert(row-1,column-1) = value;
				nnz++;
				}
			}
		}	
	M.makeCompressed();
	return M;
	}

//get last line of a file
string getLastLine(ifstream& inStream)
	{
    string xLine;
    while (inStream >> ws && getline(inStream, xLine)) // skip empty lines
    	{
        ;
		}
		
    return xLine;
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
//fourier transform type functions

//h the matrix from dl[7]
mat hFn(const unsigned int & xN, const double & xa)
	{
	mat xh(xN,xN);	xh = Eigen::MatrixXd::Zero(xN,xN);
	double diag = 1.0 + 2.0/pow(xa,2.0); //diagonal terms
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
