/* ---------------------------------------------------------------------------------------------
program to load sphaleron and sphaleronNegEig and try to calculate the zero energy instanton
in euclidean space using the newton method and an initial guess based on the solution to the
linearised equations of motion
---------------------------------------------------------------------------------------------*/
/* ---------------------------------------------------------------------------------------------
initial guess
	sphaleron plus small amount of zero mode, constant in time
---------------------------------------------------------------------------------------------*/
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
#include <ctype.h>
#include <cstring> //for memcpy
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "sphaleron_pf.h"

using namespace std;

int main(int argc, char ** argv)
{
/* ---------------------------------------------------------------------------------------------
load vectors
---------------------------------------------------------------------------------------------*/
vec sphaleron = loadSimpleVectorColumn("data/sphaleron.dat",1);
vec negEig = loadSimpleVector("data/sphaleronEigVec.dat");
if (sphaleron.size()!=negEig.size())
	{
	printf("error: sphaleron.size() = %i, negEig.size() = %i\n\n",(int)sphaleron.size(),(int)negEig.size());
	return 0;
	}

/* ---------------------------------------------------------------------------------------------
main parameters
---------------------------------------------------------------------------------------------*/
unsigned int 	N=200, Nt = 200;
double 			r0 = 1.0e-16, r1 = 10.0, t0 = 0.0, t1 = 10.0;
double			amp;
if (argc>1)
	{
	amp = atof(argv[1]);
	}
else
	{
	amp = -1.0e-2;
	}

/* ---------------------------------------------------------------------------------------------
defining main vectors
---------------------------------------------------------------------------------------------*/
vec phi(N*Nt); phi = Eigen::VectorXd::Zero(N*Nt);
vec phi0full = sphaleron + amp*negEig;
vec phi0 = interpolate1d(phi0full,sphaleron.size(),N);

/* ---------------------------------------------------------------------------------------------
constructing intial guess for phi
---------------------------------------------------------------------------------------------*/
for (unsigned int k=0; k<Nt; k++)
	{
	for (unsigned int j=0; j<N; j++)
		{
		unsigned int m = k + j*Nt;
		phi[m] = phi0[j];
		}
	}
	
/* ---------------------------------------------------------------------------------------------
printing initial guess
---------------------------------------------------------------------------------------------*/
string filename = "data/instanton00.dat";
unsigned int N_print = 300, Nt_print = 300; //these two should be equal for input to pi.cc
vec tVec(Nt_print*N_print), rVec(Nt_print*N_print), phiToPrint;
double dtPrint = (t1-t0)/(Nt_print-1.0);
double dxPrint = (r1-r0)/(N_print-1.0);
for (unsigned int k=0;k<Nt_print;k++)
	{
	for (unsigned int j=0; j<N_print; j++)
		{
		unsigned int m = k + j*Nt_print;
		tVec(m) = t0 + k*dtPrint;
		rVec(m) = r0 + j*dxPrint;
		}
	}
phiToPrint = interpolate(phi,Nt,N,Nt_print,N_print);
phiToPrint = reverseTime(phiToPrint,Nt_print,N_print);
printThreeVectors(filename,tVec,rVec,phiToPrint);
gp(filename,"vector.gp");

return 0;
}
