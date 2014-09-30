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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
//loading parameters
aqStruct aq; //struct to hold user responses

//define the numbers P and c
unsigned int P; //number of times multiplied
double c;//small number in exponent

ifstream fin;
fin.open("inputs", ios::in);
string line;
int firstLine = 0;
while(getline(fin,line))
	{
	if(line[0] == '#')
		{
		continue;
		}
	if (firstLine==0)
		{
		istringstream ss1(line);
		ss1 >> N >> Na >> Nb >> Nc >> dE >> Tb >> theta;
		firstLine++;
		}
	else if (firstLine==1)
		{
		istringstream ss2(line);
		ss2 >> aq.inputChoice >> aq.fileNo >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
		ss2 >> alpha >> open;
		firstLine++;
		}
	else
		{
		istringstream ss3(line);
		ss3 >> P >> c;
		}
	}
fin.close();
inP = aq.inputChoice; //just because I write this a lot
//derived quantities
NT = Na + Nb + Nc;
epsilon = dE;
R = 2.0/3.0/epsilon;
alpha *= R;
if (inP.compare("p") == 0)
	{
	L = 3.0*R;
	if (Tb<R)
		{
		angle = asin(Tb/R);
		double Ltemp = 1.5*(1.5*Tb*tan(angle));
		if (Ltemp<L) //making sure to use the smaller of the two possible Ls
			{
			L=Ltemp;
			}
		}
	else
		{
		cout << "Tb>R, cannot define angle" << endl;
		}
	}
else if (inP.compare("b") == 0)
	{
	Tb = 1.5*R;
	L = 3.0*R;
	}
a = L/(N-1.0);
b = Tb/(Nb-1.0);
Ta = b*Na;
Tc = b*Nc;
Gamma = exp(-theta);


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
string loadFile = "./data/DDS0" + to_string(aq.fileNo) + ".dat";

spMat M = loadSpmat(loadFile,to_reserve);

//define a random vector r
vec r = Eigen::VectorXd::Random(matSize);

//define the matrix A = (1-c*M)
spMat A(matSize,matSize);
A.setIdentity();
M = -c*M;
A += M;
M = -M/c;

//apply A to r N times in a loop - scaling in each interation
for (unsigned int j=0;j<P;j++)
	{
	r = A*r;
	double scale = r.maxCoeff();
	r /= scale;
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
temp = eigenvalue*r - temp;
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
fprintf(scalarFile,"%16g%16g%16g\n",check,error,eigenvalue);
fclose(scalarFile);

//printing eigenvector to file
string vecFile = "./data/eigVec.dat";
printVectorB(vecFile,r);

return 0;
}
