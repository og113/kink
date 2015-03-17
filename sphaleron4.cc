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
vec sphaleronFull = loadSimpleVector("data/stable/sphaleron.dat");
vec negEigFull = loadSimpleVector("data/stable/sphaleronEigVec.dat");
if (sphaleronFull.size()!=negEigFull.size()) {
	printf("error: sphaleron.size() = %i, negEig.size() = %i\n\n",(int)sphaleronFull.size(),(int)negEigFull.size());
	return 1;
}

/* ---------------------------------------------------------------------------------------------
main parameters
---------------------------------------------------------------------------------------------*/
unsigned int 	N = 1e3, Nt = 1e3;
double 			r0 = 1.0e-16, r1 = 5.0, t0 = 0.0, t1 = 0.80;
double			dr, dt;
double			amp = -1.0e-2;

/* ---------------------------------------------------------------------------------------------
getting inputs
---------------------------------------------------------------------------------------------*/
if (argc>2) {
	for (int j=0; j<(int)(argc/2); j++) {
		string temp1 = argv[2*j+1];
		string temp2 = argv[2*j+2];
		if (temp1[0]=='-') temp1 = temp1.substr(1);
		if (temp1.compare("amp")==0) amp = stringToNumber<double>(temp2);
		else if (temp1.compare("t1")==0) t1 = stringToNumber<double>(temp2);
		else if (temp1.compare("r1")==0) r1 = stringToNumber<double>(temp2);
		else if (temp1.compare("r0")==0) r0 = stringToNumber<double>(temp2);
		else if (temp1.compare("N")==0) N = stringToNumber<unsigned int>(temp2);
		else if (temp1.compare("Nt")==0) Nt = stringToNumber<unsigned int>(temp2);
		else {
			cerr << "input " << temp1 << " not understood" << endl;
			return 1;
		}
	}
}
else if (argc==2) amp = atof(argv[1]);
else cerr << "inputs not understood" << endl;

dr = (r1-r0)/(double)(N-1.0), dt = (t1-t0)/(double)(Nt-1.0);

printf("%8s%8s%8s%8s%8s%8s%8s\n","N","Nt","r1","r0","t1","t0","amp");
printf("%8i%8i%8g%8g%8g%8g%8g\n",N,Nt,r1,r0,t1,t0,amp);
printf("\n");

/* ---------------------------------------------------------------------------------------------
defining main vectors
---------------------------------------------------------------------------------------------*/
vec phi(N*Nt); phi = Eigen::VectorXd::Zero(N*Nt);
vec sphaleron = interpolate1d(sphaleronFull,sphaleronFull.size(),N);
vec negEig = interpolate1d(negEigFull,sphaleronFull.size(),N);

/* ---------------------------------------------------------------------------------------------
normalising
---------------------------------------------------------------------------------------------*/
if (negEig[0]<0) negEig *= -1.0;
double normSphaleron, normNegEig;
normSphaleron = sphaleron.norm();
normNegEig = negEig.norm();
negEig *= normSphaleron/normNegEig;

/* ---------------------------------------------------------------------------------------------
constructing intial guess for phi
---------------------------------------------------------------------------------------------*/
for (unsigned int k=0; k<Nt; k++) {
	for (unsigned int j=0; j<N; j++) {
		unsigned int 	m = k + j*Nt;
		double 			r = r0 + dr*j, t = t0 + dt*k;
		phi[m] = r*(sphaleron[j] + amp*cos(3.91*t)*negEig[j]);
	}
}
	
/* ---------------------------------------------------------------------------------------------
printing initial guess
---------------------------------------------------------------------------------------------*/
string filename = "data/instanton00.dat";
unsigned int N_print = 300, Nt_print = 300;
vec tVec(Nt_print*N_print), rVec(Nt_print*N_print), phiToPrint;
double dtPrint = (t1-t0)/(Nt_print-1.0);
double dxPrint = (r1-r0)/(N_print-1.0);
for (unsigned int k=0;k<Nt_print;k++) {
	for (unsigned int j=0; j<N_print; j++) {
		unsigned int m = k + j*Nt_print;
		tVec(m) = t0 + k*dtPrint;
		rVec(m) = r0 + j*dxPrint;
	}
}
phiToPrint = interpolate(phi,Nt,N,Nt_print,N_print);
phiToPrint = reverseTime(phiToPrint,Nt_print,N_print);
printThreeVectors(filename,tVec,rVec,phiToPrint);
//gp(filename,"repi.gp");

return 0;
}
