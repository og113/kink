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
#include <limits> //for DBL_MIN etc
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
//defining the time to label output
string timeNumber = currentDateTime();
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//solving 1d boundary value problem via shooting method
paramsVoid = {};

gsl_odeiv2_system sys = {func, jac, 4, &paramsVoid};

double F = 1.0, dF;
double aim = 0.0;
double closeness = 1.0e-8;
double r0 = 1.0e-16, r1 = 16.0;
const unsigned int N = 1e3;
double dr = r1-r0;
dr /= (double)N;
vec y0Vec(N+1), y2Vec(N+1);
unsigned int runsCount = 0;

double Y0 = 4.337;//initial guess
//cout << "initial y0: ";
//cin >> Y0;

printf("%16s%16s%16s%16s%16s%16s%16s%16s\n","run","N","y(r1)","yMin","F-aim","Y0Old","Y0New","-F/dF");

while (absolute(F-aim)>closeness)
	{
	runsCount++;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_yp_new (&sys,  gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

	double r = r0, ri;
	double y0[4] = { Y0, 0.0, 1.0, 0.0};
	double y[4];
	memcpy(y, y0, sizeof(y0));
	double yMin = absolute(y[0]-aim);
	y0Vec[0] = y[0];
	y2Vec[0] = y[2];
	int status;
	unsigned int i, iMin = 0;
	
	for (i = 1; i <= N; i++)
		{
		ri =(double)i;
		ri *= dr;
		ri += r0;
		status = gsl_odeiv2_driver_apply (d, &r, ri, y);
		if (status != GSL_SUCCESS)
			{
			printf ("error, return value=%d\n", status);
			printf ("i = %3i, r = %3g\n",i,r);
			break;
			}
		if ((y[0]-aim)<yMin && (y[0]-aim)>0.0)
			{
			yMin = y[0]-aim;
			iMin = i;
			}
		y0Vec[i] = y[0];
		y2Vec[i] = y[2];
		//printf ("%.5e %.5e %.5e\n", r, y[0], y[1]);
		if ((y[0]-aim)<0.0)
			{
			iMin = i;
			if ((y[0]-aim)<-0.2) { break; }
			}
		}
	if (status != GSL_SUCCESS){break;}
		
	F = y0Vec[iMin]-aim; //as final boundary condition is y(r1)=0.0;
	dF = y2Vec[iMin];
	printf("%16i%16i%16g%16g%16g%16.12g",runsCount,i,y[0],yMin,F-aim,Y0);
	if (absolute(dF)>2.0e-16)
		{
		Y0 += -F/dF;
		printf("%16.12g%16g\n",Y0,-F/dF);
		}
	gsl_odeiv2_driver_free (d);
	if (i==(N+1))
		{
		F = y[0];
		}
	}

//printing solution
string filename = "./data/sphaleron.dat", picname = "./pics/sphaleron.png";
simplePrintVector(filename,y0Vec);
printf("\nSolution printed: %20s %20s\n\n",filename.c_str(),picname.c_str());
gpSimple(filename,picname);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//finding E
double E, Eerror;
p  = {Y0};
gsl_function E_integrand;
E_integrand.function = &EIntegrand;
E_integrand.params = &p;
gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
gsl_integration_qag(&E_integrand, r0, r1, 1.0e-10, 1.0e-9, 1e4, 4, w, &E, &Eerror);
gsl_integration_workspace_free(w);
if (Eerror>1.0e-8) { cout << "E error = " << Eerror << endl;}
else { cout << "E = " << E << endl << endl;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//computing and printing linear fluctuation operator as sparse matrix
//in order to calulate negative eigenvalue and eigenvector
if (false)
	{
	spMat D1(N+1,N+1), D2(N+1,N+1);
	D1.setZero();
	D2.setZero();
	double r = r0;
	for (unsigned int j=0; j<(N+1); j++)
		{
		if (j==0)
			{
			D1.insert(j,j) = 1.0/pow(dr,2.0) + (1.0 - 3.0*pow(y0Vec[j],2.0))/2.0;
			D1.insert(j,j+1) = -1.0/pow(dr,2.0);
			D2.insert(j,j) = 1.0/pow(dr,2.0) + 1.0/dr/r + 1.0 - 3.0*pow(y0Vec[j],2.0);
			D2.insert(j,j+1) = -1.0/pow(dr,2.0) - 1.0/dr/r;
			/*D1.insert(j,j) = -1.0; //derivative of phi is zero at r=0
			D1.insert(j,j+1) = 1.0;
			D2.insert(j,j) = -1.0;
			D2.insert(j,j+1) = 1.0;*/
			}
		else if (j==N)
			{
			D1.insert(j,j) = 1.0/pow(dr,2.0) - 2.0*(1.0-dr/r)/r/dr + (1.0 - 3.0*pow(y0Vec[j],2.0))/2.0;
			D1.insert(j,j-1) = -1.0/pow(dr,2.0) + 2.0*(1.0-dr/r)/r/dr;
			D2.insert(j,j) = 1.0/pow(dr,2.0) - 1.0/dr/r + 1.0 - 3.0*pow(y0Vec[j],2.0);
			D2.insert(j,j-1) = -1.0/pow(dr,2.0) + 1.0/dr/r;
			/*D1.insert(j,j) = 1.0; //phi goes to zero as r->infty
			D2.insert(j,j) = 1.0;*/
			}
		else
			{
			D1.insert(j,j) = 2.0/pow(dr,2.0) - 2.0*(1.0-dr/r)/r/dr + 1.0 - 3.0*pow(y0Vec[j],2.0);
			D1.insert(j,j+1) = -1.0/pow(dr,2.0);
			D1.insert(j,j-1) = -1.0/pow(dr,2.0) + 2.0*(1.0-dr/r)/r/dr;
			D2.insert(j,j) = 2.0/pow(dr,2.0) + 1.0 - 3.0*pow(y0Vec[j],2.0);
			D2.insert(j,j+1) = -1.0/pow(dr,2.0) - 1.0/dr/r;
			D2.insert(j,j-1) = -1.0/pow(dr,2.0) + 1.0/dr/r;
			}
		r += dr;
		}
	D1.makeCompressed();
	D2.makeCompressed();
	printSpmat("./data/D1.dat",D1);
	printSpmat("./data/D2.dat",D2);
	printf("Matrices printed: ./data/D1.dat ./data/D2.dat\n\n");
	}
printf("From Matlab: D1 gives omega^2_- = -15.31,\n");
printf("             D2 gives omega^2_- = -15.34\n\n");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//assembling material to calculate occupation of modes

mat omega(N+1,N+1);
mat Eomega(N+1,N+1);
if (false)
	{
	//h_ij matrix
	mat h(N+1,N+1);
	for (unsigned int j=0;j<(N+1);j++)
		{
		if (j==0)
			{
			h(j,j) = 1.0 + 2.0/pow(dr,2.0);
			h(j,j+1) = -pow(2.0,0.5)/pow(dr,2.0);
			}
		else if (j==N)
			{
			h(j,j) = 1.0 + 2.0/pow(dr,2.0);
			h(j,j-1) = -pow(2.0,0.5)/pow(dr,2.0);
			}
		else
			{
			h(j,j) = 1.0 + 2.0/pow(dr,2.0);
			if ((j+1)==N)	{	h(j,j+1) = -pow(2.0,0.5)/pow(dr,2.0);}
			else			{	h(j,j+1) = -1.0/pow(dr,2.0);}
			if ((j-1)==0)	{	h(j,j-1) = -pow(2.0,0.5)/pow(dr,2.0);}
			else			{	h(j,j-1) = -1.0/pow(dr,2.0);}
			}
		}
	omega = Eigen::MatrixXd::Zero(N+1,N+1);
	Eomega = Eigen::MatrixXd::Zero(N+1,N+1);
	vec eigenValues(N+1);
	mat eigenVectors(N+1,N+1); //eigenvectors correspond to columns of this matrix
	Eigen::SelfAdjointEigenSolver<mat> eigensolver(h);
	if (eigensolver.info() != Eigen::Success)
		{
		cout << "h eigensolver failed" << endl;
		}
	else
		{
		eigenValues = eigensolver.eigenvalues();
		eigenVectors = eigensolver.eigenvectors(); //automatically normalised to have unit norm
		}
	//simplePrintVector("data/sphaleronhEigs.dat",eigenValues);
	
	for (unsigned int j=0; j<(N+1); j++)
		{
		for (unsigned int k=0; k<(N+1); k++)
			{
			double rj = r0 + j*dr, rk = r0 + k*dr;
			double dj = pow(4.0*pi,0.5)*rj*pow(dr,0.5);
			double dk = pow(4.0*pi,0.5)*rk*pow(dr,0.5);
			for (unsigned int l=0; l<(N+1); l++)
				{
				omega(j,k) += dj*dk*pow(eigenValues(l),0.5)*eigenVectors(j,l)*eigenVectors(k,l);
				Eomega(j,k) += dj*dk*eigenValues(l)*eigenVectors(j,l)*eigenVectors(k,l);
				}
			}
		}
	printMat("data/sphaleronOmega.dat",omega);
	printMat("data/sphaleronEomega.dat",Eomega);
	printf("omega printed: data/sphaleronOmega.dat\n");
	printf("Eomega printed: data/sphaleronOmega.dat\n\n");
	}
else
	{
	omega = loadMat("data/sphaleronOmega.dat",N+1,N+1);
	Eomega = loadMat("data/sphaleronEomega.dat",N+1,N+1);
	printf("Loaded omega and Eomega from file\n\n");
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//propagating sphaleron plus a small amount of negative mode forward in time
//in order to calculate E_lin and N_lin

double sigma = 1.0; //sigma=1 for minkowsian; sigma=-1 for euclidean
unsigned int Nt = N*3;
double T = 0.8*(r1-r0), amp = -1.0e-2;
if ((T-4.0)>1.1*(r1-r0)) { cout << "R is too small compared to T. R = " << r1-r0 << ", T = " << T << endl;}
//cout << "amp of negative mode: ";
//cin >> amp;
double dt = T/Nt; //equals Dt
if (dt>0.5*dr)
	{
	cout << "dt too large. dt = "<< dt << ", dr = " << dr << endl;
	return 0;
	}
vector<double> phi((Nt+1)*N);

vec eigVec;

//getting eigVec from file
eigVec = loadSimpleVector("data/sphaleronEigVec.dat"); //normalized to 1

paramsForce = {r0, r1, N, sigma};

//setting up ode solver for time evolution
gsl_odeiv2_system sysTime = {funcPDE, jacPDE, 4, &paramsForce};

vector<double> yVec((N+1)*(Nt+1));

gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_yp_new (&sysTime,  gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

double t = 0.0, ti;
double y0[2*N];
for (unsigned int j=0; j<N; j++)
	{
	y0[2*j] = y0Vec[j] + amp*eigVec[j];
	y0[2*j+1] = 0.0;
	phi[j*Nt] = y0[j];
	}
double y[2*N];
memcpy(y, y0, sizeof(y0));
int status;

for (unsigned int i = 1; i < Nt; i++)
	{
	ti =(double)i;
	ti *= dt;
	status = gsl_odeiv2_driver_apply (d, &t, ti, y);
	if (status != GSL_SUCCESS)
		{
		printf ("error, return value=%d\n", status);
		printf ("i = %3i, t = %3g\n",i,t);
		break;
		}
	for (unsigned int j=0; j<N; j++)
		{
		phi[i+j*Nt] = y0[2*j];
		}
	}
for (unsigned int j=0; j<N; j++)
	{
	phi[Nt+j*Nt] = 0.0;
	}
gsl_odeiv2_driver_free (d);

//printing solution
filename = "./data/sphaleron2Evo.dat";
picname = "./pics/sphaleron2Evo.png";
simplePrintVector(filename,phi);
printf("\n%16s%20s%20s\n\n","Solution printed: ",filename.c_str(),picname.c_str());
gpSimple(filename,picname);

return 0;
}
