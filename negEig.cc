//program to find the smallest eigenvector and eigenvalue for a large sparse matrix
//used in conjuction with pi.cc
//loads from file a triplet (i,j,v) defining the sparse matrix
//prints to files the eigenvector,  eigenvalue and error
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
#include <gsl/gsl_poly.h>
#include "gnuplot_i.hpp"
#include "pf.h"
#include "files.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
string timeNumber;
cout << "input timeNumber: ";
cin >> timeNumber;

//load parameters
//getting variables and user inputs from inputs
ifstream fin;
fin.open("inputs");
if (fin.is_open())
	{
	string line;
	unsigned int lineNumber = 0;
	while(getline(fin,line))
		{
		if(line[0] == '#')
			{
			continue;
			}
		istringstream ss(line);
		if (lineNumber==0)
			{
			ss >> N >> Na >> Nb >> Nc >> dE >> LoR >> Tb >> theta;
			lineNumber++;
			}
		else if (lineNumber==1)
			{
			ss >> aq.inputChoice >> aq.inputFile >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
			lineNumber++;
			}
		else if(lineNumber==2)
			{
			ss >> alpha >> open >> amp >> pot >> A >> reg;
			lineNumber++;
			}
		else if(lineNumber==3)
			{
			double temp;
			ss >> temp >> temp >> negEigDone;
			lineNumber++;
			}
		}
	}
else
	{
	cout << "unable to open inputs" << endl;
	}
fin.close();
inP = aq.inputChoice; //just because I write this a lot

//derived quantities
NT = Na + Nb + Nc;
epsilon = dE;
R = 2.0/3.0/epsilon;
L = LoR*R;
if (inP.compare("p") == 0)
	{
	if (Tb<R)
		{
		angle = asin(Tb/R);
		double Ltemp = 1.5*(1.5*Tb*tan(angle));
		if (Ltemp<L) //making sure to use the smaller of the two possible Ls
			{
			L=Ltemp;
			}
		}
	}
else if (inP.compare("b") == 0)
	{
	Tb = 1.5*R;
	}
a = L/(N-1.0);
b = Tb/(Nb-1.0);

//define the numbers P and c
unsigned int P; //number of times multiplied
double c;//small number in exponent

ifstream inputs;
inputs.open("inputs", ios::in);
string lastLine = getLastLine(inputs);
istringstream ss(lastLine);
ss >> P >> c;

//load DDS from file into matrix M
unsigned int zeroModes = 1;
unsigned int matSize = 2*N*Nb+zeroModes;
Eigen::VectorXi to_reserve(matSize);
to_reserve = Eigen::VectorXi::Constant(matSize,11);
to_reserve(0) = 3; //these need to be changed when boundary conditions need to be more compicated
to_reserve(1) = 3;
to_reserve(2*N*Nb-2) = 3;
to_reserve(2*N*Nb-1) = 3;
to_reserve(2*N*Nb) = N;
unsigned int fileNumber = 0;
string loadFile = "./data/" + timeNumber + "DDSb" + "_" + numberToString<unsigned int>(fileNumber) + ".dat";

spMat M = loadSpmat(loadFile,to_reserve);

//define random vectors r, s, t
vec r = Eigen::VectorXd::Random(matSize);
vec s = Eigen::VectorXd::Random(matSize);
vec t = Eigen::VectorXd::Random(matSize);

//define the matrix A = (1-c*M)
spMat A(matSize,matSize);
A.setIdentity();
M = -c*M;
A += M;
M = -M/c;

//apply A to r N times in a loop - scaling in each interation
unsigned int P1 = floor(P/10);
for (unsigned int j=0;j<P1;j++)
	{
	r = A*r;
	s = A*s;
	t = A*t;
	double scaleR = r.maxCoeff();
	double scaleS = s.maxCoeff();
	double scaleT = t.maxCoeff();
	r /= scaleR;
	s /= scaleS;
	t /= scaleT;
	}
r += s + t;
r /= 3.0;

for (unsigned int j=P1; j<P;j++)
	{
	r = A*r;
	double scaleR = r.maxCoeff();
	r /= scaleR;
	}

//normalising r to get approximation to eigenvector v
double norm = r.dot(r);
norm = pow(norm,0.5);
r /= norm;

//evaluating (v,Mv) to get approximation to eigenvalue lambda
double eigenvalue;
vec temp = M*r;
eigenvalue = r.dot(temp);

//evaluating error as |lambda*v-M*v|
double error;
temp += -eigenvalue*r;
error = temp.dot(temp);
error = pow(error,0.5);

//checking that epsilon*lambda*N<1 - is this really what we need
double check = absolute(c*eigenvalue*N);
if (check>1)
	{
	//cout << "check>1" << endl;//not sure about the relevance of this
	}

//printing error, and eigenvalue to file
FILE * scalarFile;
scalarFile = fopen("./data/eigValue.dat","a");
fprintf(scalarFile,"%16i%16g%16g%16g%16g\n",P,c,check,error,eigenvalue);
fclose(scalarFile);

//printing eigenvector to file
string vecFile = "./data/eigVec.dat";
printVectorB(vecFile,r);

return 0;
}
