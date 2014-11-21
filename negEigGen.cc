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
//taking in timeNumber
string filename;
cout << endl;
cout << "filename: ";
cin >> filename;

//starting clock
clock_t time;
time = clock();

//define the numbers P and c
unsigned int P; //number of times multiplied
double c;//small number in exponent

ifstream inputs;
inputs.open("inputs", ios::in);
string lastLine = getLastLine(inputs);
istringstream ss(lastLine);
ss >> P >> c;
cout << "P: " << P << endl;
cout << "c: " << c << endl;
cout << endl;

//load matrix from file into matrix M
spMat M = cleverLoadSpmat(filename);
unsigned int matSize = M.rows();

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

//checking that epsilon*lambda*M<1 - is this really what we need
double check = absolute(c*eigenvalue*P);
//if (check>1)
	//{
	//cout << "check>1" << endl;//not sure about the relevance of this
	//}
	
//stopping clock
time = clock() - time;
double realtime = time/1000000.0;
	
//printing time, error, and eigenvalue to screen
cout << left << setw(20) << "time" << setw(20) << "check<1" << setw(20) << "error" << setw(20) << "eigenvalue" << endl;
cout << left << setw(20) << realtime << setw(20) << check << setw(20) << error << setw(20) << eigenvalue << endl;
cout << endl;

//printing error, and eigenvalue to file
FILE * scalarFile;
scalarFile = fopen("./data/eigValue.dat","a");
fprintf(scalarFile,"%16i%16g%16g%16g%16g\n",P,c,check,error,eigenvalue);
fclose(scalarFile);

//printing eigenvector to file
//string vecFile = "./data/eigVec.dat";
//printVectorB(vecFile,r);

return 0;
}
