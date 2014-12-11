/* first of all just writing something to solve the ode as an initial value plot
only once that works will i add the newton raphson loop to do the shooting
to do:
- define gsl functions for first order odes
- define jacobian functions
- define boundary values
- follow gsl example for implementation

next, have to implement shooting by applying newton's method to F, the difference between the final boundary value and the one shot at
- define function F
- write loop with convergence test
- write newton's method
- print tests for convergence on each iteration

then,
- calculate energy of solution
- find negative mode
- calculate number of particles
	- to do this need to evolve sphaleron + small amount of negative mode in time
	- then model N(k) and E(K) and N_lin and E_lin with time
	- this requires taking something like the fourier transform of the field configuration
*/
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
const unsigned int N = 2e2;
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
string filename = "./data/" + timeNumber + "sphaleron.dat", picname = "./pics/sphaleron.png";
simplePrintVector(filename,y0Vec);
printf("\n%16s%20s%20s\n\n","Solution printed: ",filename.c_str(),picname.c_str());
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

unsigned int Nt = N*3;
double T = 0.8*(r1-r0), amp = -1.0e-2, sigma = 1.0; //set sigma=-1 for euclidean evolution
if ((T-4.0)>1.1*(r1-r0)) { cout << "R is too small compared to T. R = " << r1-r0 << ", T = " << T << endl;}
//cout << "amp of negative mode: ";
//cin >> amp;
double dt = T/Nt; //equals Dt
if (dt>0.5*dr)
	{
	cout << "dt too large. dt = "<< dt << ", dr = " << dr << endl;
	return 0;
	}
vec phi((Nt+1)*(N+1)), vel((Nt+1)*(N+1)), acc((Nt+1)*(N+1)), eigVec, linNum(Nt+1), linErg(Nt+1), linErgField(Nt+1);
vec nonLinErg(Nt+1), erg(Nt+1);
linErg = Eigen::VectorXd::Zero(Nt+1);
linErgField = Eigen::VectorXd::Zero(Nt+1);
linNum = Eigen::VectorXd::Zero(Nt+1);
nonLinErg = Eigen::VectorXd::Zero(Nt+1);
erg = Eigen::VectorXd::Zero(Nt+1);

//getting eigVec from file
eigVec = loadSimpleVector("data/sphaleronEigVec.dat"); //normalized to 1

if (eigVec.size()!=(N+1)){ printf("eigVec.size() = %i. Redo Matlab calculation of eigVec.\n\n",(int)(eigVec.size())); return 0;}

//initial condition 1)
//defining phi on initial time slice
for (unsigned int j=0;j<(N+1);j++)
	{
	unsigned int l = j*(Nt+1);
	phi(l) = y0Vec[j] + amp*eigVec(j);
	}
	
//initialising linErg and linNum	
for (unsigned int x=0;x<(N+1);x++)
	{
	unsigned int j = x*(Nt+1);
	double r = r0 + x*dr, sigma = 1.0;
	if (x==0 || x==N) {sigma = 0.5;}
	nonLinErg(0) += 4.0*pi*pow(r,2.0)*sigma*0.25*pow(phi(j),4.0)*dr;
	linErgField(0) += 4.0*pi*pow(r,2.0)*sigma*0.5*pow(phi(j),2.0)*dr;
	if (x<N)
		{
		linErgField(0) += 4.0*pi*r*(r+dr)*0.5*pow(phi(j+(Nt+1))-phi(j),2.0)/dr;
		}
	for (unsigned int k=0;k<(N+1);k++)
		{
		unsigned int l = k*(Nt+1);
		linErg(0) += Eomega(x,k)*phi(l)*phi(j);
		linNum(0) += omega(x,k)*phi(l)*phi(j);	
		}
	}

//initial condition 2)
//intitialize velocity, dphi/dt(x,0) = 0;
for (unsigned int j=0; j<(N+1); j++)
	{
	unsigned int l = j*(Nt+1);
    vel(l) = 0.0;
	}

//boundary condition 1)
//phi=0.0 at r=L
for (unsigned int j=0;j<(Nt+1);j++)
	{
	unsigned int l = N*(Nt+1) + j;
	phi(l) = 0.0;
	}

//boundary condition 2)
//initialize acc using phi and expression from equation of motion
/*the unusual d^2phi/dr^2 term and the absence of the first derivative term
are due to the boundary condition 2) which, to second order, is phi(t,-dr) = phi(t,dr)*/
acc(0) = 2.0*(phi(Nt+1) - phi(0))/pow(dr,2.0) - phi(0) + pow(phi(0),3.0);
acc(0) *= sigma*0.5; //as initial time slice, generated from taylor expansion and equation of motion
for (unsigned int j=1; j<N; j++)
	{
	unsigned int l = j*(Nt+1);
	double r = r0 + dr*j;
    acc(l) = (phi(l+(Nt+1)) + phi(l-(Nt+1)) - 2.0*phi(l))/pow(dr,2.0) + (phi(l+(Nt+1))-phi(l-(Nt+1)))/r/dr - phi(l) + pow(phi(l),3.0);
    acc(l) *= sigma*0.5;
	}

//A7. run loop
for (unsigned int u=1; u<(Nt+1); u++)
	{
    for (unsigned int x=0; x<N; x++) //don't loop over last x position as fixed by boundary condition 1)
    	{
        unsigned int m = u+x*(Nt+1);
        vel(m) = vel(m-1) + dt*acc(m-1);
        phi(m) = phi(m-1) + dt*vel(m);
    	}
    acc(u) = 2.0*(phi(u+Nt+1) - phi(u))/pow(dr,2.0) - phi(u) + pow(phi(u),3.0);
    acc(u) *= sigma;
    linErgField(u-1) += 4.0*pi*dr*pow(r0,2.0)*0.5*pow(phi(u)-phi(u-1),2.0)/pow(dt,2.0);
    linErgField(u-1) += 4.0*pi*dr*pow(r1,2.0)*0.5*pow(phi(u+N*(Nt+1))-phi(u+N*(Nt+1)-1),2.0)/pow(dt,2.0);
    linErgField(u) += 4.0*pi*( r0*(r0+dr)*0.5*pow(phi(u+(Nt+1))-phi(u),2.0)/dr + dr*pow(r0,2.0)*0.5*pow(phi(u),2.0) );
    linErgField(u) += 4.0*pi*dr*pow(r1,2.0)*0.5*pow(phi(u+N*(Nt+1)),2.0);
    nonLinErg(u) += 4.0*pi*dr*pow(r0,2.0)*0.25*pow(phi(u),4.0);
    nonLinErg(u) += 4.0*pi*dr*pow(r1,2.0)*0.25*pow(phi(u+N*(Nt+1)),4.0);
    for (unsigned int x=1; x<N; x++)
    	{
        unsigned int m = u+x*(Nt+1);
        double r = r0 + x*dr;
        acc(m) = (phi(m+(Nt+1)) + phi(m-(Nt+1)) - 2.0*phi(m))/pow(dr,2.0) + (phi(m+(Nt+1))-phi(m-(Nt+1)))/r/dr - phi(m) + pow(phi(m),3.0);
        acc(m) *= sigma;
        linErgField(u-1) +=  4.0*pi*pow(r,2.0)*dr*0.5*pow(phi(m)-phi(m-1),2.0)/pow(dt,2.0);
		linErgField(u) += 4.0*pi*(r*(r+dr)*0.5*pow(phi(m+(Nt+1))-phi(m),2.0)/dr + pow(r,2.0)*0.5*dr*pow(phi(m),2.0));
		nonLinErg(u) += 4.0*pi*pow(r,2.0)*0.25*pow(phi(m),4.0)*dr;
        for (unsigned int k=0;k<(N+1);k++)
			{
			unsigned int j = u+k*(Nt+1);
			linErg(u) += Eomega(x,k)*phi(m)*phi(j);
			linNum(u) += omega(x,k)*phi(m)*phi(j);
			}
    	}
	}
	
for (unsigned int k=0; k<(Nt+1); k++)
	{
	erg(k) = linErgField(k) + nonLinErg(k);
	}
	
double linErgContm = 0.0, linNumContm = 0.0;
for (unsigned int k=1; k<(N+1); k++)
	{
	double momtm = k*pi/(r1-r0);
	double freqSqrd = 1.0+pow(momtm,2.0);
	double Asqrd, integral1 = 0.0, integral2 = 0.0;
	for (unsigned int l=0; l<(N+1); l++)
		{
		double r = r0 + l*dr;
		unsigned int m = (Nt-1) + l*(Nt+1);
		integral1 += dr*r*phi(m)*pow(2.0/(r1-r0),0.5)*sin(momtm*r);
		integral2 += dr*r*(phi(m+1)-phi(m))*pow(2.0/(r1-r0),0.5)*sin(momtm*r)/dt;
		}
	Asqrd = pow(integral1,2.0) + pow(integral2,2.0)/freqSqrd;
	linErgContm += 2.0*pi*Asqrd*freqSqrd;
	linNumContm += 2.0*pi*Asqrd*pow(freqSqrd,0.5);
	}
	
unsigned int N_print = 100, Nt_print = 100;
vec tVec(Nt_print*N_print), rVec(Nt_print*N_print);
double dtPrint = T/(Nt_print-1.0);
double dxPrint = (r1-r0)/(N_print-1.0);
for (unsigned int t=0;t<Nt_print;t++)
	{
	for (unsigned int r=0; r<N_print; r++)
		{
		unsigned int j= t + r*Nt_print;
		tVec(j) = t*dtPrint;
		rVec(j) = r0 + r*dxPrint;
		}
	}
	
vec phiToPrint;
phiToPrint = interpolate(phi,Nt+1,N+1,Nt_print,N_print);
simplePrintVector("data/sphaleronLinNum.dat",linNum);
simplePrintVector("data/sphaleronLinErg.dat",linErg);
simplePrintVector("data/sphaleronLinErgField.dat",linErgField);
simplePrintVector("data/sphaleronNonLinErg.dat",nonLinErg);
simplePrintVector("data/sphaleronErg.dat",erg);
gpSimple("data/sphaleronLinNum.dat","pics/sphaleronLinNum.png");
gpSimple("data/sphaleronLinErg.dat","pics/sphaleronLinErg.png");
gpSimple("data/sphaleronLinErgField.dat","pics/sphaleronLinErgField.png");
gpSimple("data/sphaleronNonLinErg.dat","pics/sphaleronNonLinErg.png");
gpSimple("data/sphaleronErg.dat","pics/sphaleronErg.png");
string evoPrint = "data/" + timeNumber + "sphaleronEvo.dat";
printThreeVectors(evoPrint,tVec,rVec,phiToPrint);
gp(evoPrint,"sphaleron.gp");
printf("Time evolution printed: %39s pics/sphaleronEvolution.png\n",evoPrint.c_str());
printf("                        data/sphaleronLinNum.dat pics/sphaleronLinNum.png\n");
printf("                        data/sphaleronLinErg.dat pics/sphaleronLinErg.png\n");
printf("                        data/sphaleronLinErgField.dat pics/sphaleronLinErgField.png\n");
printf("                        data/sphaleronNonLinErg.dat pics/sphaleronNonLinErg.png\n");
printf("                        data/sphaleronErg.dat pics/sphaleronErg.png\n\n");
printf("erg(Nt-1) = %8.4f\n",linErgField(Nt-1));
printf("linErgField(Nt-1) = %8.4f\n",linErgField(Nt-1));
printf("nonLinErg(Nt-1) = %8.4f\n",nonLinErg(Nt-1));
printf("linNum(Nt-1) = %8.4f\n",linNum(Nt-1));
printf("linErg(Nt-1) = %8.4f\n",linErg(Nt-1));
printf("linNumContm(Nt-1) = %8.4f\n",linNumContm);
printf("linErgContm(Nt-1) = %8.4f\n\n",linErgContm);

return 0;
}
