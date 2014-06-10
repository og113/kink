#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

complex<double> i(0.0,1.0);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//primary parameters
unsigned int N = 64; //number of points in each spatial dimension
unsigned int Na = (int)(1.2*N);
unsigned int Nb = (int)(1.0*N);
unsigned int Nc = 2;
double R = 40.0; //size of bubble
double mass = 4.0; 
double lambda = 0.1;
double Tb = 1.2*R/2;
double angle = arcsin(Lb/R);
double L = 2*(1.5*Lb*tan(angle));

//derived quantities
unsigned int Nt = Na + Nb + Nc;
double a = L/(N-1.0); //step sizes in each spatial dimension
double b = Tb/(Nb-1.0); //step sizes in time
double Ta = b*(Na-1.0);
double Tc = b*(Nc-1.0);
double v =  mass*pow(lambda,-0.5); //vacuum phi
double X = mass*R; //the large thin-wall parameter, the ratio of the size of the bubble to the size of the wall
double epsilon = 2.0*pow(mass,3)/lambda/R/3.0;; //energy difference

//determining number of runs
double closeness = pow(10,-3.0);
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

//print string with short spaces
void coutShort (const vector<string>& labels)
	{
	cout << left;
	for (unsigned l=0; l<labels.size(); l++)
		{
		cout << setw(8) << labels[l];
		}
	cout << endl;
	}
	
//print string with long spaces
void coutLong (const vector<string>& labels)
	{
	cout << left;
	for (unsigned l=0; l<labels.size(); l++)
		{
		cout << setw(16) << labels[l];
		}
	cout << endl;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//askQuestions and changeParameters

//struct to hold answers to questions
struct aq{
	unsigned int fileNo;
	double maxTheta;
	unsigned int totalLoops;
	string printChoice;
	unsigned int printRun;
};


