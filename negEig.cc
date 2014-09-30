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
//load DDS from file into matrix M
unsigned int Nt, zeroModes, fileNo;

ifstream fin;
fin.open("negEigInputs", ios::in);
string line;
while(getline(fin,line))
	{
	if(line[0] == '#')
		{
		continue;
		}
	istringstream ss(line);
	ss >> N >> Nt >> zeroModes >> fileNo;
	}
fin.close();

unsigned int matSize = 2*N*Nt+zeroModes;
unsigned int nnzGuess = (Nt-2)*(2*N*11) + 6 + 6 + N*zeroModes; //change Nt to suit, and change 6 to a maximum of 2*N for more complicated boundaries
string loadFile = "./data/DDS0" + to_string(fileNo) + ".dat";

spMat M = loadSpmat(loadFile,nnzGuess);

//define a random vector r
vec r = Eigen::VectorXd::Random(matSize);

//define the numbers P and c
unsigned int P = 100; //number of times multiplied
double c = 0.01;

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
double check = c*eigenvalue*N;
if (check<1)
	{
	cout << "check<1" << endl;
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
