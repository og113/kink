/* something to compile results to get periodic instanton for pot3
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

int main(int argc, char** argv)
{
//defining the time to label output
string timeNumber = currentDateTime();
string timeNumberIn = "150112114306";
int direction = 1;
double sigma = 1.0;
if (argc == 2) timeNumberIn = argv[1];
else if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id.compare("-tn")==0) timeNumberIn = argv[2*j+2];
		else if (id.compare("-d")==0) direction = atoi(argv[2*j+2]);
		else if (id.compare("-s")==0) sigma = atof(argv[2*j+2]);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//main parameters

double r0 = 1.0e-16, r1 = 10.0;
const unsigned int N = 1e3;
double dr = r1-r0;
dr /= (double)N;
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//propagating sphaleron plus a small amount of negative mode forward in time
//in order to calculate E_lin and N_lin

unsigned int Nt = N*4;
double T = (double)direction*0.8*(r1-r0); //set sigma=-1 for euclidean evolution
if ((abs(T)-4.0)>1.1*(r1-r0)) {
	cout << "R is too small compared to T. R = " << r1-r0 << ", T = " << T << endl;
	return 0;
}
//cout << "amp of negative mode: ";
//cin >> amp;
double dt = T/Nt; //equals Dt
if (abs(dt)>0.5*dr) {
	cout << "dt too large. dt = "<< dt << ", dr = " << dr << endl;
	return 0;
}
vec phi((Nt+1)*(N+1)), vel((Nt+1)*(N+1)), acc((Nt+1)*(N+1)), linErgField(Nt+1);
vec nonLinErg(Nt+1), erg(Nt+1);
linErgField = Eigen::VectorXd::Zero(Nt+1);
nonLinErg = Eigen::VectorXd::Zero(Nt+1);
erg = Eigen::VectorXd::Zero(Nt+1);

//getting initial from file
string filename = "data/" + timeNumberIn + "pip_0.dat";
unsigned int fileLength = countLines(filename);
vec initial;
		
if (fileLength==((N+1)*(Nt+1)))
	{
	initial = loadSimpleVectorColumn(filename,3);
	}
else if (fileLength % 2) //if its odd
	{
	unsigned int Nin, Ntin;
	//cout << "interpolating input, filelength = " << fileLength << " , Cp.size() = " << N*Nb+1 << endl;
	string inputsF = "./data/" + timeNumberIn + "inputsPi_0";
	ifstream fin;
	fin.open(inputsF.c_str());
	if (fin.is_open())
		{
		string line, temp;
		while(getline(fin,line))
			{
			if(line[0] == '#') continue;
			istringstream ss(line);
			ss >> Nin >> temp >> Ntin;
			break;
			}
		}
	else cout << "unable to open " << inputsF << endl;
	fin.close();
	vec temp_initial = loadSimpleVectorColumn(filename,3);
	temp_initial.conservativeResize(Ntin*Nin);
	initial = interpolate(temp_initial,Ntin,Nin,Nt+1,N+1);
	}

//initial condition 1)
//defining phi on initial time slice
for (unsigned int j=0;j<(N+1);j++)
	{
	unsigned int l = j*(Nt+1);
	unsigned int m;
	(direction == 1) ? m = Nt + l : m = l;
	double r = r0 + j*dr;
	phi(l) = initial(m)/r;
	}
	
//initialising linErg and linNum	
for (unsigned int x=0;x<(N+1);x++)
	{
	unsigned int j = x*(Nt+1);
	double r = r0 + x*dr, eta;
	(x==0 || x==N) ? eta = 0.5 : eta = 1.0;
	nonLinErg(0) += 4.0*pi*pow(r,2.0)*eta*0.25*pow(phi(j),4.0)*dr;
	linErgField(0) += 4.0*pi*pow(r,2.0)*eta*0.5*pow(phi(j),2.0)*dr;
	if (x<N) linErgField(0) += 4.0*pi*r*(r+dr)*0.5*pow(phi(j+(Nt+1))-phi(j),2.0)/dr;
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

string evoPrint = "data/" + timeNumber + "pi3.dat";
printThreeVectors(evoPrint,tVec,rVec,phiToPrint);
gp(evoPrint,"pi3.gp");
printf("\ndirection = %i, sigma = %g\n",direction,sigma);
printf("Input:                  %39s\n",filename.c_str());
printf("Time evolution printed: %39s pics/pi3.png\n",evoPrint.c_str());
printf("erg(Nt-1) = %8.4f\n",linErgField(Nt-1));
printf("linErgField(Nt-1) = %8.4f\n",linErgField(Nt-1));
printf("nonLinErg(Nt-1) = %8.4f\n",nonLinErg(Nt-1));
printf("linNumContm(Nt-1) = %8.4f\n",linNumContm);
printf("linErgContm(Nt-1) = %8.4f\n\n",linErgContm);

return 0;
}
