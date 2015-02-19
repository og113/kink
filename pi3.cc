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

/* ---------------------------------------------------------------------------------------------
main parameters
---------------------------------------------------------------------------------------------*/

double r0 = 1.0e-16, r1 = 10.0;
unsigned int N = 300;
double dr;

double Ta, Tc;
unsigned int Na, Nc, Nt;
double dt;

int direction = 1; // direction of time evolution
double sigma = 1.0; //set sigma=-1 for euclidean evolution
bool testTunnel = false, testLinear = false, changeParams = false, approxOmega = false;
double closenessLin = 1.0e-2, Tlin = 0.0;
double closenessEdge = 1.0e-3, Redge = 0.0;
double closenessMom = 1.9e-2, momTest;

/* ---------------------------------------------------------------------------------------------
user inputs
and labels for input and output
---------------------------------------------------------------------------------------------*/

string timeNumber;
string timeNumberIn = "150112114306";
string loopIn = "0";
if (argc == 2) timeNumberIn = argv[1];
else if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("tn")==0) timeNumberIn = argv[2*j+2];
		else if (id.compare("r1")==0) timeNumberIn = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("loop")==0 || id.compare("l")==0) loopIn = argv[2*j+2];
		else if (id.compare("test")==0 || id.compare("tunnel")==0 || id.compare("tt")==0) testTunnel = (bool)atoi(argv[2*j+2]);
		else if (id.compare("linearization")==0 || id.compare("lin")==0) testLinear = (bool)atoi(argv[2*j+2]);
		else if (id.compare("closeness")==0 || id.compare("close")==0) closenessLin = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("changeParams")==0) changeParams = (bool)atoi(argv[2*j+2]);
		else if (id.compare("approxOmega")==0) approxOmega = (bool)atoi(argv[2*j+2]);
		else if (id.compare("N")==0) N = atoi(argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}
else if (argc == 1) {
	cerr << "must provide input to pi3.cc" << endl;
	return 1;
}
else {
	cerr << "must provide an even number of inputs in format '-name value':" << endl;
	for (int j=1; j<argc; j++) cerr << argv[j] << " ";
	cerr << endl;
	return 1;
}

timeNumber = timeNumberIn;
if (testTunnel) {
	cout << "pi3 testing if tunnelled" << endl;
}
if (testLinear) {
	testTunnel = false;
	cout << "pi3 testing linearization" << endl;
}

/* ---------------------------------------------------------------------------------------------
getting parameters from specific inputs
---------------------------------------------------------------------------------------------*/

unsigned int Nin, Nain, Nbin, Ncin;
double Tbin, dtin, drin;
string inputsF = "./data/" + timeNumberIn + "inputsPi_" + loopIn;
ifstream fin;
fin.open(inputsF.c_str());
if (fin.is_open())
	{
	string line, temp;
	double LoR;
	while(getline(fin,line))
		{
		if(line[0] == '#') continue;
		if(line.empty()) continue;
		istringstream ss(line);
		ss >> Nin >> Nain >> Nbin >> Ncin >> temp >> LoR >> Tbin;
		r1 = 10.0*LoR;
		break;
		}
	}
else cout << "unable to open " << inputsF << endl;
fin.close();

dr = r1-r0;
dr /= (double)N;
dtin = Tbin/(Nbin-1.0);
drin = (r1-r0)/(Nin-1.0);
Ta = dtin*Nain;
Tc = dtin*Ncin;
dt = dr*0.2;
Na = (unsigned int)(Ta/dt);
Nc = (unsigned int)(Tc/dt);

if (abs(Ta)>1.1*(r1-r0) && !testTunnel && !testLinear) {
	cerr << "R is too small compared to Ta. R = " << r1-r0 << ", Ta = " << Ta << endl;
	return 1;
}
if (abs(Tc)>1.1*(r1-r0) && !testTunnel && !testLinear) {
	cerr << "R is too small compared to Tc. R = " << r1-r0 << ", Tc = " << Tc << endl;
	return 1;
}
if (abs(dt)>0.5*dr) {
	cerr << "dt too large. dt = " << dt << ", dr = " << dr << endl;
	return 1;
}
if (abs(Ta)<2.5 && !testTunnel && !testLinear) {
	cerr << "Ta too small. Ta = " << Ta << endl;
	return 1;
}

/* ---------------------------------------------------------------------------------------------
loading phi on BC
---------------------------------------------------------------------------------------------*/
string filename = "data/" + timeNumberIn + "pip_" + loopIn + ".dat";
unsigned int fileLength = countLines(filename);
vec phiBC;
	
if (fileLength==Nin*Nbin) phiBC = loadSimpleVectorColumn(filename,3);
else if (fileLength==(Nin*Nbin+1))
	{
	phiBC = loadSimpleVectorColumn(filename,3);
	phiBC.conservativeResize(Nbin*Nin);
	}
else {
	cerr << filename << " cannot be read." << endl;
	cerr << "File length not correct: " << fileLength << " != " << Nbin*Nin+1 << endl;
}


/* ---------------------------------------------------------------------------------------------
deterimining omega matrices for fourier transforms in spatial direction
---------------------------------------------------------------------------------------------*/
mat omega_m1(N+1,N+1); 	
mat omega_0(N+1,N+1);
mat omega_1(N+1,N+1);
mat omega_2(N+1,N+1);
vec eigenValues(N+1);
mat eigenVectors(N+1,N+1); //eigenvectors correspond to columns of this matrix
if (testLinear) {
	omega_m1 = Eigen::MatrixXd::Zero(N+1,N+1);
	omega_0 = Eigen::MatrixXd::Zero(N+1,N+1);
	omega_1 = Eigen::MatrixXd::Zero(N+1,N+1);
	omega_2 = Eigen::MatrixXd::Zero(N+1,N+1);
	if(!approxOmega) {
		mat h(N+1,N+1);
		h = hFn(N+1,dr,1.0);
		Eigen::SelfAdjointEigenSolver<mat> eigensolver(h);
		if (eigensolver.info() != Eigen::Success) cerr << "h eigensolver failed" << endl;
		else {
			eigenValues = eigensolver.eigenvalues();
			eigenVectors = eigensolver.eigenvectors(); //automatically normalised to have unit norm
		}
	}
	else {
		double normalisation = sqrt(2.0/(double)N);
		for (unsigned int l=0; l<(N+1); l++) {
			eigenValues(l) = 1.0+pow(2.0*sin(pi*l/(double)N/2.0)/dr,2.0);
			for (unsigned int m=0; m<(N+1); m++) eigenVectors(l,m) = normalisation*sin(pi*l*m/(double)N);
		}
	}
	double djdk;	
	for (unsigned int j=0; j<(N+1); j++) {
		for (unsigned int k=0; k<(N+1); k++) {
			for (unsigned int l=0; l<(N+1); l++) {
				djdk = 4.0*pi*dr;
				if ((j==0 || j==N) && !approxOmega) djdk/=sqrt(2.0);
				if ((k==0 || k==N) && !approxOmega) djdk/=sqrt(2.0);
				omega_m1(j,k) += djdk*pow(eigenValues(l),-0.5)*eigenVectors(j,l)*eigenVectors(k,l);
				omega_0(j,k) += djdk*eigenVectors(j,l)*eigenVectors(k,l);
				omega_1(j,k) += djdk*pow(eigenValues(l),0.5)*eigenVectors(j,l)*eigenVectors(k,l);
				omega_2(j,k) += djdk*eigenValues(l)*eigenVectors(j,l)*eigenVectors(k,l);
			}
		}
	}
}

/* ---------------------------------------------------------------------------------------------
propagating euclidean solution forwards and/or backwards in time
---------------------------------------------------------------------------------------------*/

vec phiA, phiC, linearizationA;
double linErgContm, linNumContm, nonLinErgA, linErgFieldA, ergA, linErgA, linNumA;
double kineticT, kineticS, massTerm;
double Tinf = 0.0;

uint j=0;
while(j<2) {

	if (testLinear) {
		Nt = 10*N;
		j=1;
		direction = -1;
	}
	else if (testTunnel) {
		Nt = 10*N;
		j=1;
		direction = 1;
	}
	else if (j==0) {
		Nt = Nc;
		direction = 1;
	}
	else if (j==1) {
		Nt = Na;
		direction = -1;
	}
	
	vec phi((Nt+1)*(N+1)), vel((Nt+1)*(N+1)), acc((Nt+1)*(N+1));
	vec nonLinErg(Nt+1), linErgField(Nt+1), erg(Nt+1);
	linErgField = Eigen::VectorXd::Zero(Nt+1);
	nonLinErg = Eigen::VectorXd::Zero(Nt+1);
	erg = Eigen::VectorXd::Zero(Nt+1);
	linErgContm = 0.0, linNumContm = 0.0;
	linErgA = 0.0, linNumA = 0.0;
	kineticT = 0.0, kineticS = 0.0, massTerm = 0.0;
	vec initial;
	initial = interpolate(phiBC,Nbin,Nin,Nt+1,N+1);

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
		eta = ((x==0 || x==N) ? 0.5 : 1.0);
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
		for (unsigned int x=0; x<N; x++) { //don't loop over last x position as fixed by boundary condition 1)
		    unsigned int m = u+x*(Nt+1);
		    vel(m) = vel(m-1) + dt*acc(m-1);
		    phi(m) = phi(m-1) + dt*vel(m);
		    if (testTunnel) {
		    	double testInf = phi(m);
		    	if (!isfinite(testInf) && abs(Tinf)<1.0e-16) Tinf = u*dt;
		    }
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
	
	for (unsigned int k=0; k<(N+1); k++) {
		unsigned int u = Nt + k*(Nt+1);
		double r = r0 + k*dr;
		kineticT += 4.0*pi*dr*pow(r,2.0)*0.5*pow(phi(u)-phi(u-1),2.0)/pow(dt,2.0);
		if (k<N) kineticS += 4.0*pi*r*(r+dr)*0.5*pow(phi(u+(Nt+1))-phi(u),2.0)/dr;
		massTerm += 4.0*pi*dr*pow(r,2.0)*0.5*pow(phi(u-1),2.0);
	}
	
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
	
	if (j==0) {
		phiC = phi;
	}
	else if (j==1) {
		phiA = phi;
		ergA = erg(Nt-1);
		nonLinErgA = nonLinErg(Nt-1);
		linErgFieldA = linErgField(Nt-1);
		if (testLinear) {
			linearizationA = Eigen::VectorXd::Zero(Nt+1);
			bool nonLin = true;
			for (unsigned int k=0; k<(Nt+1); k++) {
				linearizationA(k) = absDiff(erg(k),linErgField(k));
				if (linearizationA(k)>closenessLin) nonLin = true;
				if (linearizationA(k)<closenessLin && nonLin) {
					Tlin = k*dt;
					nonLin = false;
				}
			}
			for (unsigned int j=0;j<(N+1);j++){
				unsigned int n = (N-j)*(Nt+1)+(unsigned int)(Tlin/dt);
				double rj = r0+(N-j)*dr;
				if (abs(phi(n))>closenessEdge && abs(Redge)<1.0e-16)
		        	Redge = rj;
				for (unsigned int k=0;k<(N+1);k++){
					unsigned int l = j*(Nt+1)+Nt-1;
					unsigned int m = k*(Nt+1)+Nt-1;
					double r = r0 + j*dr;
					double s = r0 + k*dr;
					// the would be imaginary parts of the following expression cancel as omega and Eomega are symmetric
					// so wlog we have dropped them
						linErgA += 0.5*r*s*\
							(omega_2(j,k)*phi(l)*phi(m)+omega_0(j,k)*(phi(l+1)-phi(l))*(phi(m+1)-phi(m))/pow(dt,2.0));
						linNumA += 0.5*r*s*\
							(omega_1(j,k)*phi(l)*phi(m)+omega_m1(j,k)*(phi(l+1)-phi(l))*(phi(m+1)-phi(m))/pow(dt,2.0));
				}
			}
		}
	}
	j++;
		
} // end of while j<2 loop


if (testLinear) {
	vec tVec = Eigen::VectorXd::Zero(Nt+1);
	for (unsigned int t=0;t<(Nt+1);t++){
		tVec(t) = t*dt;
		}
	string linearizationFile = "data/" + timeNumber + "linearization_" + loopIn + ".dat";
	simplePrintVector(linearizationFile,tVec);
	simpleAppendVector(linearizationFile,linearizationA);
	momTest = linErgA*dr/linNumA/pi;
	printf("linearization printed:  %39s\n",linearizationFile.c_str());
	printf("Tlin(%6.4f)          = %6.4f\n",closenessLin,Tlin);
	printf("Redge(%6.4f)         = %6.4f\n",closenessEdge,Redge);
	printf("momTest               = %6.4f\n",momTest);
	if (changeParams) {
		//string inputsOut = "data/" + timeNumber + "inputsPi3_" + loopIn;
		string inputsOut = "inputs";
		int Nmom = (int)(linErgA*Redge/pi/linNumA/closenessMom) + 1;
		int Nbmom = (int)(Nmom*4*Tbin/Redge);
		int Nalin = (int)(Nbmom*Tlin/Tbin);
		printf("L changed to :        %6.4g\n",Redge);
		printf("N changed to :        %6i\n",Nmom);
		printf("Na changed to:        %6i\n",Nalin); 
		printf("Nb changed to:        %6i\n\n",Nbmom); 
		changeInputs("data/temp1","LoR",numberToString<double>(Redge/10.0),inputsF);
		changeInputs("data/temp2","N",numberToString<int>(Nmom),"data/temp1");
		changeInputs("data/temp1","Na",numberToString<int>(Nalin),"data/temp2");
		changeInputs("data/temp2","Nb",numberToString<int>(Nbmom),"data/temp1");
		copyFile("data/temp2",inputsOut);
	}
}
else if (!testTunnel) {	
	vec tVec((Nain+Nbin+Ncin)*Nin), rVec((Nain+Nbin+Ncin)*Nin);
	for (unsigned int t=0;t<(Nain+Nbin+Ncin);t++)
		{
		for (unsigned int r=0; r<Nin; r++)
			{
			unsigned int j= t + r*(Nain+Nbin+Ncin);
			tVec(j) = t*dtin;
			rVec(j) = r0 + r*drin;
			}
		}
	
	// constructing input to main
	vec mainIn((Nain+Nbin+Ncin)*Nin);
	vec phiAOut, phiCOut;
	phiAOut = interpolate(phiA,Na+1,N+1,Nain+1,Nin);
	phiCOut = interpolate(phiC,Nc+1,N+1,Ncin+1,Nin);
	for (unsigned int j=0;j<(Nain+Nbin+Ncin);j++)
		{
		for (unsigned int k=0; k<Nin; k++)
			{
			unsigned int l = j+k*(Nain+Nbin+Ncin), m;
			double r = r0 + k*drin;
			if (j<Nain)
				{
				m = (Nain-j)+k*(Nain+1);
				mainIn(l) = phiAOut(m)*r;
				}
			else if (j<(Nain+Nbin))
				{
				m = (j-Nain)+k*Nbin;
				mainIn(l) = phiBC(m);
				}
			else
				{
		        m = j - Nain - Nbin + 1 + k*(Ncin+1);
				mainIn(l) = phiCOut(m)*r;
				}
			}
		}
	string mainInFile = "data/" + timeNumber + "tpip_" + loopIn + ".dat";
	printThreeVectors(mainInFile,tVec,rVec,mainIn);
	gp(mainInFile,"pi3.gp");

	printf("%8s%8s%8s%8s%8s%8s\n","N","Na","Nb","Nc","L","Tb");
	printf("%8i%8i%8i%8i%8g%8g\n",Nin,Nain,Nbin,Ncin,r1-r0,Tbin);
	printf("\n");
	
	printf("Input:                  %39s\n",filename.c_str());
	printf("tpip printed:                %39s pics/pi3.png\n",mainInFile.c_str());
}
else {
	printf("Input:              %39s\n",filename.c_str());
	printf("Tinf              = %8.4f\n",Tinf);
	}
				printf("kineticT          = %8.4f\n",kineticT);
				printf("kineticS          = %8.4f\n",kineticS);
				printf("massTerm          = %8.4f\n",massTerm);
				printf("erg(0)            = %8.4f\n",ergA);
				printf("nonLinErgA(0)     = %8.4f\n",nonLinErgA);
				printf("linErgFieldA(0)   = %8.4f\n",linErgFieldA);
				printf("linErgContmA(0)   = %8.4f\n",linErgContm);
if (testLinear) printf("linErgA           = %8.4f\n",linErgA);
				printf("linNumContmA(0)   = %8.4f\n",linNumContm);
if (testLinear) printf("linNumA           = %8.4f\n\n",linNumA);

double finalTest = linErgFieldA;
if (testTunnel) {
	if ( !isfinite(finalTest) ) return 0;
	else return 1;
}
else return 0;
}
