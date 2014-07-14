#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

typedef unsigned long int lint;
typedef complex <double> comp;
typedef vector <unsigned int> intVec;
typedef vector <comp> cVec;
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
double R = 10.0; //size of bubble
double mass = 3.0; 
double lambda = 0.1;
double Tb = 1.2*R/2;
double angle = arcsin(Tb/R); %not a primary parameter, just used to make L
double L = 2*(1.5*Tb*tan(angle));

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
	return Xintcoord;	
	}
	
//gives values of coordinates in whole spacetime
complex<double> coord(const unsigned& int locNum,const int& direction)
	{
	complex<double> Xcoord;
	if (direction==0)
		{
		unsigned int t = intCoords(locNum,0,Nt);
		if ( t < Na)
			{
			Xcoord = b*(t-Na) + i*Tb;
			}
		else if (intCoords(locNum,0) < (Na+Nb))
			{
			Xcoord = i*Tb - i*b*(t-Na);
			}
		else
			{
			Xcoord = b*(t-Na-Nb);
			}
		}
	if (direction==1)
		{
		unsigned int x = intCoords(locNum,1,Nt);
		Xcoord = -L/2.0 + x*a;
		}
	return Xcoord;
	}

//gives values of coordinates on section AB
complex<double> coordA(const unsigned int& locNum,const int& direction)
	{
	complex<double> XcoordA;
	if (direction==0)
		{
		unsigned int t = intCoords(locNum,0,Na);
		XcoordA = b*(t-Na) + i*Tb;
		}
	if (direction==1)
		{
		unsigned int x = intCoords(locNum,1,Na);
		XcoordA = -L/2.0 + x*a;
		}
	return XcoordA;
	}

//gives values of coordinates on section BC
complex<double> coordB(const unsigned int& locNum,const int& direction)
	{
	complex<double> XcoordB;
	if (direction==0)
		{
		unsigned int t = intCoords(locNum,0,Nb);
		XcoordB = i*(Tb - b*t);
		}
	if (direction==1)
		{
		unsigned int x = intCoords(locNum,1,Nb);
		XcoordB = -L/2.0 + x*a;
		}
	return XcoordB;
	}
	
//gives values of coordinates on section CD
complex<double> coordC(const unsigned int& locNum,const int& direction)
	{
	complex<double> XcoordC;
	if (direction==0)
		{
		unsigned int t = intCoords(locNum,0,Nc);
		XcoordC = b*t;
		}
	if (direction==1)
		{
		unsigned int x = intCoords(locNum,1,Nc);
		XcoordC = -L/2.0 + x*a;
		}
	return XcoordC;
	}
	
long int neigh(const lint& locNum, const unsigned int& direction, const signed int& sign, const unsigned int& xNt) //periodic in space but not time, degree refers to the number of neighbours, 1 is for just positive neighbours, 2 is for both
	{
	long int neighLocation = -1; //this is the result if there are no neighbours for the given values of the argument
	unsigned int c = intCoords(locNum,direction,xNt);
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
	fprintf("%8s","N","Na","Nb","Nc","L","Tb","R","mass","lambda","theta");
	fprintf("\n");
	fprintf("%8g",N,Na,Nb,Nc,L,Tb,R,mass,lambda,theta);
	fprintf("\n");
	}
	
//print action and its constituents to the terminal
void printAction ( const comp& Kinetic, const comp& potL, const comp& potE)
	{
	action = Kinetic + potL + potE;
	fprintf("%16s","kinetic","potL","potE","action");
	fprintf("\n");
	fprintf("%16g",Kinetic,potL,potE,action);
	fprintf("\n");
	}
	
//print vector to file
void printVector (const string& printFile, vec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
			for (unsigned long int j=0; j<Eucdim; j++)
				{
				F << left;
				for (int r=0; r<dim; r++)
					{
					F << setw(15) << re(coord(j,r)) << setw(15) << im(coord(j,r));
					}
				F << setw(15) << vecToPrint(2*j) << setw(15) << vecToPrint(2*j+1)  << endl;		
				}
			F.close();
	}
	
//print sparse matrix to file
void printSpmat (const string & printFile, spMat spmatToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F << left;
	for (unsigned int l=0; l<spmatToPrint.outerSize(); l++)
		{
		for (SparseMatrix<double>::InnerIterator it(spmatToPrint,l); it; it++)
			{
			F << setw(15) << it.row() << setw(15) << it.col() << setw(15) << it.value() << endl;
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
	string loopParameter;
	double minValue;
	double maxValue;
	string printChoice;
	unsigned int printRun;
	};

//asks initial questions to get user inputs
void askQuestions (aqStruct * aqx )
	{
	cout << "number of loops: ";
	cin >> aqx.totalLoops;
	cout << endl;
	if aqx.totalLoops !=1
		{
		cout << "loop (main)parameter (N,Na,Nb,Nc,L,Tb,R,mass,lambda) ";
		cin >> aqx.loopParameter;
		cout << end;
		cout << "min value): ";
		cin >> aqx.minValue;
		cout << end;
		cout << "max value: ";
		cin >> aqx.maxValue;
		cout << end;
		}
	cout << "print early (m,v,p,a,n)? ";
	cin >> aqx.printChoice;
	cout << end;
	cout << "which run to print? ";
	cin >> aqx.printRun	
	}

//changes parameters according to user inputs (integer)
void changeInt (const string & parameterLabel, const int & newParameter)
	{
	if ( paramLabel.compare("N")==0)
		{
		Na = (int)newParameter*Na/N;
		Nb = (int)newParameter*Nb/N;
		Nc = (int)newParameter*Nc/N;
		Nt = Na + Nb + Nc;
		N = newParameter;
		a = L/(N-1);
		b = Tb/(Nb-1);
		Ta = b*(Na-1.0);
		Tc = b*(Nc-1.0);
		}
	elseif ( paramLabel.compare("Na")==0)
		{
		Na = newParameter;
		Nt = Na + Nb + Nc;
		Ta = b*(Na-1.0);
		}
	elseif ( paramLabel.compare("Nb")==0)
		{
		Nb = newParameter;
		Nt = Na + Nb + Nc;
		b = Tb/(Nb-1.0);
		Ta = b*(Na-1.0);
		Tc = b*(Nc-1.0);
		}
	elseif ( paramLabel.compare("Nc")==0)
		{
		Nc = newParameter;
		Nt = Nc + Nb + Nc;
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
	if ( paramLabel.compare("L")==0)
		{
		L = newParameter;
		a = L/(N-1);
		}
	else if ( paramLabel.compare("Tb")==0)
		{
		Tb = newParameter;
		b = Tb/(Nb-1.0);
		Ta = b*(Na-1.0);
		Tc = b*(Nc-1.0);
		angle = arcsin(Tb/R);
		L = 2*(1.5*Tb*tan(angle));
		}
	else if ( paramLabel.compare("R")==0)
		{
		R = newParameter;
		X = mass*R;
		epsilon = 2.0*pow(mass,3)/lambda/R/3.0;
		angle = arcsin(Tb/R);
		L = 2*(1.5*Tb*tan(angle));
		}
	else if ( paramLabel.compare("mass")==0)
		{
		mass = newParameter;
		v =  mass*pow(lambda,-0.5);
		X = mass*R;
		epsilon = 2.0*pow(mass,3)/lambda/R/3.0;
		}
	else if ( paramLabel.compare("lambda")==0)
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
	vec real(2*tDim);
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
void convergence (const int& runsCount, const double& actionTest, const double& vectorTest, const clock_t& Clock, clock_t* Wait, char* printWait, char* printChoice, int* printRun, const comp& Kinetic, const comp& potL, const comp& potE, bool boolWait)
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
			double Action = Kinetic+potL+potE;
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
	
