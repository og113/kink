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

using namespace std;

typedef unsigned long int lint;
typedef complex<double> comp;
typedef vector<unsigned int> intVec;
typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;

int main()
{
//load DDS from file into matrix M

//define a random vector r

//define the numbers N and epsilon

//define the matrix A = (1-epsilon*M)

//apply A to r N times in a loop

//normalising r to get approximation to eigenvector v

//evaluating (v,Mv) to get approximation to eigenvalue lambda

//evaluating error as |lambda*V-M*v|

//checking that epsilon*lambda*N<1

//printing error, and eigenvalue to file

//printing eigenvector to file
return 0;
}
