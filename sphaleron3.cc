/* ---------------------------------------------------------------------------------------------
program to load sphaleron and sphaleronNegEig and try to calculate the zero energy instanton
in euclidean space using the newton method and an initial guess based on the solution to the
linearised equations of motion
---------------------------------------------------------------------------------------------*/
/* ---------------------------------------------------------------------------------------------
initial guess
	initial condition expanded in solutions of the linear pde
	i.e. spherical bessel functions in r multiplied by decaying exponentials in t
	the intial guess then satisfies both boundary conditions and the initial condition
	but not the final condition, though it gets exponentially close for large t1-t0
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

int main()
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
if (negEig[0]<0) negEig *= -1.0;

/* ---------------------------------------------------------------------------------------------
main parameters
---------------------------------------------------------------------------------------------*/
unsigned int 	N=200, Nt = 200;
double 			r0 = 1.0e-16, r1 = 10.0, t0 = 0.0, t1 = 10.0;
double			dr = (r1-r0)/(double)(N-1.0), dt = (t1-t0)/(double)(Nt-1.0);
double			closeness = 1.0e-6;
double			amp = -1.0e-2;

/* ---------------------------------------------------------------------------------------------
defining main vectors
---------------------------------------------------------------------------------------------*/
vec phi(N*Nt); phi = Eigen::VectorXd::Zero(N*Nt);
vec phi0full = sphaleron + amp*negEig;
vec phi0 = interpolate1d(phi0full,sphaleron.size(),N);

/* ---------------------------------------------------------------------------------------------
finding coefficients of decomposition into spherical bessel functions, sin(kr)/kr
---------------------------------------------------------------------------------------------*/ 
unsigned int cutoff = (unsigned int)((N-1)/1);
vec coefficients(cutoff); coefficients = Eigen::VectorXd::Zero(cutoff);
for (unsigned int n=1; n<(cutoff+1); n++)
	{
	for (unsigned int j=1; j<N; j++)
		{
		double rj = r0 + j*dr;
		double waveNoj = pi*j/(double)(N);
		coefficients[n-1] += (2.0/(double)(N)) * rj * sin(waveNoj*n) * phi0[j];
		}
	}
	
/* ---------------------------------------------------------------------------------------------
testing completeness of basis
---------------------------------------------------------------------------------------------*/ 
bool completeness = true;
if (completeness)
	{
	vec testCompleteness(N); testCompleteness = Eigen::VectorXd::Zero(N);
	double completenessError = 0.0;
	for (unsigned int j=0; j<N; j++)
		{
		if (j==0)
			{
			testCompleteness[j] = phi0[j];
			}
		else
			{
			for (unsigned int k=1; k<(cutoff+1); k++)
				{
				double rj = r0 + j*dr;
				double waveNok = pi*k/(double)(N);
				testCompleteness[j] += coefficients[k-1] * ( sin(waveNok*j)/rj );
				}
			}
		testCompleteness[j] -= phi0[j];
		if (absolute(testCompleteness[j])>completenessError) { completenessError = absolute(testCompleteness[j]);}
		}
	if (completenessError>closeness)
		{
		string filename = "data/completenessError.dat", picname = "pics/completenessError.png";
		printf("\nCompleteness test maximum error: %g > %g\n",completenessError,closeness);
		simplePrintVector(filename,testCompleteness);
		gpSimple(filename,picname);
		printf("Completeness error vector printed: %20s %20s\n\n",filename.c_str(),picname.c_str());
		}
	}
	
/* ---------------------------------------------------------------------------------------------
constructing intial guess for phi
---------------------------------------------------------------------------------------------*/
for (unsigned int k=0; k<Nt; k++)
	{
	double tk = t0 + k*dt;
	for (unsigned int j=0; j<N; j++)
		{
		double rj = r0 + j*dr;
		unsigned int m = k + j*Nt;
		if (j==0)
			{
			phi[m] = phi0[j] * exp( - tk );
			}
		else
			{
			for (unsigned int l=1; l<(cutoff+1); l++)
				{
				double waveNol = pi*l/(double)(N);
				phi[m] += coefficients[l-1] * ( sin(waveNol*j)/rj ) * exp( -pow(1.0 + pow(l*pi/(r1-r0),2.0),0.5)*tk );
				}
			}
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
