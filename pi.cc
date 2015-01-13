//program to generate the periodic instantons
//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works
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
#include <limits>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include "gnuplot_i.hpp"
#include "pf.h"
#include "files.h"

using namespace std;

int main()
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//getting variables and user inputs from inputs

//defining the time to label output
bool printTimeNumber = true;
string timeNumber;
if (printTimeNumber) timeNumber = currentDateTime();

ifstream fin;
fin.open("inputs");
if (fin.is_open())
	{
	string line;
	unsigned int lineNumber = 0;
	while(getline(fin,line))
		{
		if(line[0] == '#')
			{
			continue;
			}
		istringstream ss(line);
		if (lineNumber==0)
			{
			ss >> N >> Na >> Nb >> Nc >> dE >> LoR >> Tb >> theta;
			lineNumber++;
			if (absolute(theta)>1.0e-16)
				{
				cout << "theta != 0" << endl;
				cout << "theta = " << theta << endl;
				}
			}
		else if (lineNumber==1)
			{
			ss >> aq.inputChoice >> aq.inputTimeNumber >> aq.inputLoop >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
			lineNumber++;
			}
		else if(lineNumber==2)
			{
			ss >> alpha >> open >> amp >> pot >> A >> reg;
			lineNumber++;
			}
		else if(lineNumber==3)
			{
			double temp;
			ss >> temp >> temp >> negEigDone;
			lineNumber++;
			}
		}
	}
else
	{
	cout << "unable to open inputs" << endl;
	}
fin.close();
inP = aq.inputChoice; //just because I write this a lot

string loop_choice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
string print_choice = aq.printChoice;
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//potential functions
	
	//potential functions
	if (pot[0]=='1')
		{
		neigh = &periodic;
		V = &V1c;
		dV = &dV1c;
		ddV = &ddV1c;
		Vd = &V1;
		dVd = &dV1;
		ddVd = &ddV1;
		epsilon0 = 0.0;
		epsilon = dE;
		}
	else if (pot[0]=='2')
		{
		neigh = &periodic;
		V = &V2c;
		dV = &dV2c;
		ddV = &ddV2c;
		Vd = &V2;
		dVd = &dV2;
		ddVd = &ddV2;
		epsilon0 = 0.74507774287199924;
		epsilon = 0.75;
		}
	else if (pot[0]=='3')
		{
		neigh = &spherical;
		V = &V3c;
		dV = &dV3c;
		ddV = &ddV3c;
		Vd = &V3;
		dVd = &dV3;
		ddVd = &ddV3;
		epsilon0 = 0.0;
		epsilon = 0.0;
		r0 = 1.0e-16;
		}
	else
		{
		cout << "pot option not available, pot = " << pot << endl;
		}
	paramsV  = {epsilon, A};
	paramsV0 = {epsilon0, A};
	paramsVoid = {};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//finding epsilon and root
	vector<double> minima0(2);
	
	if (pot[0]!='3')
		{
		//gsl function for dV(phi)
		gsl_function F;
		F.function = Vd;
		F.params = &paramsV;	
	
		//finding preliminary roots of dV(phi)=0
		minima[0] = brentMinimum(&F, -1.1, -3.0, 0.0);
		minima[1] = brentMinimum(&F, 1.2, 0.5, 3.0);
	
		//gsl function for V(root2)-V(root1)-dE
		struct ec_params ec_params = { A, minima[0], minima[1], dE};
		gsl_function EC;
		EC.function = &ec;
		EC.params = &ec_params;

		//evaluating epsilon, new root and dE may change slightly
		epsilonFn(&F,&EC,&dE,&epsilon,&minima);
	
		//evaluating some properties of V
		mass2 = ddVd(minima[0],&paramsV);
		cout << "mass2 = " << mass2 << endl;
		double fn_mass2 = 1.0 + epsilon0*exp(-4.0/pow(A,2.0))*(1.0/A + 1.0/pow(A,3.0) + 1.0/pow(A,5.0));
		cout << "fn_mass2 = " << fn_mass2 << endl;
	
		//finding root0 of dV0(phi)=0;
		if (pot[0]=='1')
			{
			minima0[0] = -1.0; minima0[1] = 1.0;
			}
		else if (pot[0]=='2')
			{
			gsl_function V0;
			V0.function = Vd;
			V0.params = &paramsV0;	
			minima0[0] = brentMinimum(&V0, -1.0, -3.0, 0.0);
			minima0[0] = brentMinimum(&V0, 1.2, 0.5, 3.0);
			struct ec_params ec0_params = { A, minima0[0], minima0[1], 0.0};
			gsl_function EC0;
			EC0.function = &ec;
			EC0.params = &ec0_params;
			double dE0 = 0.0;
			epsilonFn(&V0,&EC0,&dE0,&epsilon0,&minima0);
			}
	
		//finding S1
		double S1error;
		gsl_function S1_integrand;
		S1_integrand.function = &s1Integrand;
		S1_integrand.params = &paramsV0;
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
		gsl_integration_qag(&S1_integrand, minima0[0], minima0[1], 1.0e-16, 1.0e-8, 1e4, 4, w, &S1, &S1error);
		gsl_integration_workspace_free(w);
		if (S1error>1.0e-8) { cout << "S1 error = " << S1error << endl;}
		}
	else
		{
		mass2 = 1.0;
		minima[0] = 0.0; //only one minimum
		R = 1.0; alpha = 0.0; //not used
		}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//finding phi profile between minima
	
	unsigned int profileSize = Nb; //more than the minimum
	vector<double> phiProfile(profileSize);
	vector<double> rhoProfile(profileSize);
	double alphaL = alpha, alphaR = alpha;

	if (pot[0]!='3')
		{
		if (pot[0]=='2')
			{
			double phiL = minima0[1]-1.0e-2;
			double phiR = minima0[0]+1.0e-2;
			for (unsigned int j=0;j<profileSize;j++)
				{
				phiProfile[j] = phiL + (phiR-phiL)*j/(profileSize-1.0);
				}
	
			double profileError;
			gsl_function rho_integrand;
			rho_integrand.function = &rhoIntegrand;
			rho_integrand.params = &paramsV0;
			gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
			w = gsl_integration_workspace_alloc(1e4);
			for (unsigned int j=0;j<profileSize;j++)
				{
				gsl_integration_qags(&rho_integrand, phiProfile[j], 0, 1.0e-16, 1.0e-6, 1e4, w, &(rhoProfile[j]), &profileError);
				if (profileError>1.0e-5) { cout << "profile error = " << profileError << " , for j = " << j << endl;}
				}
			gsl_integration_workspace_free(w);
			alphaL = rhoProfile[0];
			alphaR = rhoProfile.back();
			}
		}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//other derived quantities
NT = Na + Nb + Nc;
Gamma = exp(-theta);
double twaction;
if (pot[0]!='3')
	{
	R = S1/dE;
	twaction = -pi*epsilon*pow(R,2)/2.0 + pi*R*S1;
	alpha *= R;
	L = LoR*R;
	}
else
	{
	twaction = 8.0*pow(pi,2.0)/3.0; //massless limit
	L = LoR*10.0;
	//no need for R or alpha
	}
vec negVec(2*N*Nb+1);
unsigned int negP;
double negc;
double negcheck;
double negerror; //should be <<1
if ((inP.compare("p") == 0 || inP.compare("f") == 0) && pot[0]!='3')
	{
	if (Tb<R)
		{
		angle = asin(Tb/R);
		double Ltemp = 1.5*(1.5*Tb*tan(angle));
		if (Ltemp<L) //making sure to use the smaller of the two possible Ls
			{
			L=Ltemp;
			}
		}
	else
		{
		//cout << "Tb>R so using negEig, need to have run pi with inP='b'" << endl;
		if (negEigDone==0)
			{
			cout << "need to run negEig and set negEigDone=1" << endl;
			return 1;
			}
		unsigned int fileLength = countLines("./data/eigVec.dat");
		if (fileLength==(N*Nb+1))
			{
			negVec = loadVector("./data/eigVec.dat",Nb,N,1);
			}
		else
			{
			cout << "eigVec not the right length" << endl;
			return 1;
			}
		ifstream eigFile;
		eigFile.open("./data/eigVal.dat", ios::in);
		string lastLine = getLastLine(eigFile);
		istringstream ss(lastLine);
		ss >> negP >> negc >> negcheck >> negerror >> negVal;
		eigFile.close();
		if (negerror>1.0)
			{
			cout << "error in negEig = " << negerror << endl;
			cout << "consider restarting with different values of P and c" << endl;
			}
		}
	}
else if ((inP.compare("b") == 0) && pot[0]!='3')
	{
	Tb = 1.5*R;
	}
else
	{
	//negVec = loadVector("data/sphaleron
	}
a = L/(N-1.0);
b = Tb/(Nb-1.0);
if (a>pow(mass2,0.5) || b>pow(mass2,0.5)) {cout << endl << "a = " << a << " , b = " << b << endl << endl;}
Ta = b*Na;
Tc = b*Nc;
double ergZero = N*a*Vd(minima[0],&paramsV);

//determining number of runs
closenessA = 1.0;
closenessS = 1.0e-5;
closenessSM = 1.0e-4;
closenessD = 1.0;
closenessC = 1.0e-16*N*NT;
closenessE = 1.0e-2;
closenessR = 1.0e-2;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//begin loop over varying parameter
//(NB. may only be one loop)
//and defining some quantities to be used later
for (unsigned int loop=0; loop<aq.totalLoops; loop++)
	{
	//giving values of varying parameters
	int intLoopParameter;
	double doubleLoopParameter;
	if (loop_choice[0] == 'N')
		{
		intLoopParameter = aq.minValue + (int)(aq.maxValue - aq.minValue)*loop/(aq.totalLoops-1);
		changeInt (loop_choice,intLoopParameter);
		}
	else if (loop_choice.compare("n")!=0)
		{
		doubleLoopParameter = aq.minValue + (aq.maxValue - aq.minValue)*loop/(aq.totalLoops-1.0);
		changeDouble (loop_choice,doubleLoopParameter);
		}
		
	//defining a time and starting the clock
	clock_t time;
	time = clock();
	
	//printing loop name and parameters
	printf("%12s%12s\n","timeNumber: ",timeNumber.c_str());
	printParameters();
	//printMoreParameters();
	
	comp action = twaction;
	double W;
	double E;
	cVec erg(NT);	

	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp action_last = action;
	unsigned int runs_count = 0;
	unsigned int min_runs = 2;
	vector<double> action_test(1);	action_test[0] = 1.0;
	vector<double> sol_test(1);		sol_test[0] = 1.0;
	vector<double> solM_test(1);	solM_test[0] = 1.0;
	vector<double> delta_test(1); 	delta_test[0] = 1.0;
	vector<double> calc_test(1); 	calc_test[0] = 1.0;
	vector<double> erg_test(1); 	erg_test[0] = 1.0;
	vector<double> reg_test(1); 	reg_test[0] = 1.0;

	//initializing phi (=p)
	vec p(2*N*Nb+1);
	p = Eigen::VectorXd::Zero(2*N*Nb+1);
	
	//defining lambda functions for regularization
	auto Vr = [&] (const comp & phi)
		{
		return -i*reg*VrFn(phi,minima[0],minima[1]);
		};
	auto dVr = [&] (const comp & phi)
		{
		return -i*reg*dVrFn(phi,minima[0],minima[1]);
		};
	auto ddVr = [&] (const comp & phi)
		{
		return -i*reg*ddVrFn(phi,minima[0],minima[1]);
		};

	//deterimining omega matrices for fourier transforms in spatial direction
	mat h(N,N);
	h = hFn(N,a,mass2);
	mat omega(N,N); 	omega = Eigen::MatrixXd::Zero(N,N);
	mat Eomega(N,N); 	Eomega = Eigen::MatrixXd::Zero(N,N);
	vec eigenValues(N);
	mat eigenVectors(N,N); //eigenvectors correspond to columns of this matrix
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
	//#pragma omp parallel for	
	for (unsigned int j=0; j<N; j++)
		{
		for (unsigned int k=0; k<N; k++)
			{
			for (unsigned int l=0; l<N; l++)
				{
				double djdk;
				if (pot[0]=='3')
					{
					double rj = r0 + j*a, rk = r0 + k*a;
					djdk = 4.0*pi*rj*rk*a;			
					}
				else
					{
					djdk=a;
					}
				if (j==0 || j==N || k==0 || k==N) {djdk/=2.0;}
				omega(j,k) += djdk*pow(eigenValues(l),0.5)*eigenVectors(j,l)*eigenVectors(k,l);
				Eomega(j,k) += djdk*eigenValues(l)*eigenVectors(j,l)*eigenVectors(k,l);
				}
			}
		}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//assigning input phi
	if (inP.compare("f")!=0 && loop==0)
		{
		if (pot[0]=='3')
			{
			vec tempPhi = loadVectorColumn("data/instanton00.dat",2);
			vec temp2Phi;
			unsigned int length = tempPhi.size();
			length = (unsigned int)(sqrt(length));
			temp2Phi = interpolate2(tempPhi,length,length,Nb,N);
			for (unsigned int j=0; j<N*Nb; j++)
				{
				p[2*j] = temp2Phi[j];
				p[2*j+1] = 0.0;
				}
			p[2*Nb*N] = 0.5; //for the zero mode
			}
		else
			{
			if (R<alpha)
				{
				cout << "R is too small. Not possible to give thinwall input. It should be more that " << alpha;
				}
			//#pragma omp parallel for
			for (unsigned int j=0; j<N*Nb; j++)
				{
				comp t = coordB(j,0);
				comp x = coordB(j,1);
				p(2*j+1) = 0.0; //imaginary parts set to zero
				if (inP.compare("b")==0 || (inP.compare("p")==0 && Tb>R))
					{
					double rho = real(sqrt(-pow(t,2.0) + pow(x,2.0))); //should be real even without real()
					if (pot[0]=='1')
						{
						if ((rho-R)<-alpha)
							{
							p(2*j) = minima[1];
							}
						else if ((rho-R)>alpha)
							{
							p(2*j) = minima[0];
							}
						else
							{
							p(2*j) = (minima[1]+minima[0])/2.0 + (minima[0]-minima[1])*tanh((rho-R)/2.0)/2.0;
							}
						}
					else if (pot[0]=='2')
						{
						if ((rho-R)<=alphaL)
							{
							p(2*j) = minima[1];
							}
						else if ((rho-R)>=alphaR)
							{
							p(2*j) = minima[0];
							}
						else
							{
							vector<double> rhoPos (profileSize,rho-R);
							for (unsigned int k=0; k<profileSize; k++)
								{
								rhoPos[k] -= rhoProfile[k];
								}
							unsigned int minLoc = smallestLoc(rhoPos);
			                p(2*j) = phiProfile[minLoc];
							}
						}
					if (inP.compare("p")==0)
						{
						p(2*j) += amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j);
						p(2*j+1) += amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j+1);
						}
					}
				else if (inP.compare("p")==0 && Tb<R)
					{
					double rho1 = real(sqrt(-pow(t,2.0) + pow(x+R*cos(angle),2.0)));
					double rho2 = real(sqrt(-pow(t,2.0) + pow(x-R*cos(angle),2.0)));
					if ((rho1-R)<-alpha && (rho2-R)<-alpha)
						{
						p(2*j) = minima[1];
						}
					else if ((rho1-R)>alpha || (rho2-R)>alpha)
						{
						p(2*j) = minima[0];
						}
					else if (real(x)>0) //note that the coord should be real
						{
						p(2*j) = (minima[1]+minima[0])/2.0 + (minima[0]-minima[1])*tanh((rho1-R)/2.0)/2.0;
						}
					else if (real(x)<0)
						{
						p(2*j) = (minima[1]+minima[0])/2.0 + (minima[0]-minima[1])*tanh((rho2-R)/2.0)/2.0;
						}
					else
						{
						p(2*j) = minima[1]; //i.e. if coordB(j,1) == 0
						}
					}
				}
			}
		p(2*N*Nb) = 0.5; //initializing Lagrange parameter for removing dp/dx zero mode
		}
	else
		{
		string loadfile;
		if (inP.compare("f")==0)
			{
			loadfile = "./data/" + aq.inputTimeNumber + "pip_" + aq.inputLoop + ".dat";
			cout << "input: " << loadfile << endl;
			inP = "p";
			}
		else
			{
			loadfile = "./data/" + timeNumber + "pi"+inP+"_" + numberToString<int>(loop-1)+".dat";
			}
		unsigned int fileLength = countLines(loadfile);
		if (fileLength==(N*Nb+1))
			{
			p = loadVector(loadfile,Nb,N,1);
			}
		else if (fileLength % 2) //if its odd
			{
			unsigned int Nin, Ntin;
			//cout << "interpolating input, filelength = " << fileLength << " , Cp.size() = " << N*Nb+1 << endl;
			string inputsF = "./data/" + aq.inputTimeNumber + "inputsPi_" + aq.inputLoop;
			ifstream fin;
			fin.open(inputsF.c_str());
			if (fin.is_open())
				{
				string line, temp;
				while(getline(fin,line))
					{
					if(line[0] == '#')
						{
						continue;
						}
					istringstream ss(line);
					ss >> Nin >> temp >> Ntin;
					break;
					}
				}
			else{cout << "unable to open " << inputsF << endl;}
			fin.close();
			vec temp_p = loadVector(loadfile,Ntin,Nin,1);
			p = interpolate(temp_p,Ntin,Nin,Nb,N);
			}
		else
			{
			cout << "vector in file: " << loadfile << " has length " << fileLength << " so cannot load into p" << endl;
			}
		}
	
	//fixing input periodic instanton to have zero time derivative at time boundaries
    if (pot[0]!='3')
    	{
    	//#pragma omp parallel for
		for (unsigned int j=0;j<N;j++)
			{
		    p(2*j*Nb) = (1.0-open)*p(2*j*Nb) + open*p(2*(j*Nb+1)); //initial time real
		    p(2*(j*Nb+1)) = p(2*j*Nb);
		    p(2*j*Nb+1) = (1.0-open)*p(2*j*Nb+1) + open*p(2*(j*Nb+1)+1); //initial time imag
		    p(2*(j*Nb+1)+1) = p(2*j*Nb+1);
		    p(2*((j+1)*Nb-1)) = open*p(2*((j+1)*Nb-2)) + (1.0-open)*p(2*((j+1)*Nb-1)); //final time real
		    p(2*((j+1)*Nb-2)) = p(2*((j+1)*Nb-1));
		    p(2*((j+1)*Nb-2)+1) = open*p(2*((j+1)*Nb-1)+1) + (1.0-open)*p(2*((j+1)*Nb-2)+1); //final time imag
		    p(2*((j+1)*Nb-1)+1) = p(2*((j+1)*Nb-2)+1);
			}
		}
	else
		{
		for (unsigned int j=0;j<Nb;j++)
			{
			unsigned int m = j + (N-1)*Nb;
		    p(2*m) = 0.0; //p=0 ar r=R
			}
		}
		
	//very early vector print
	string earlyPrintFile = "data/" + timeNumber + "piE"+inP+ "_" + numberToString<unsigned int>(loop) + "_0.dat";
	printVectorB(earlyPrintFile,p);
		
	//defining complexified vector Cp
	cVec Cp(Nb*N);
	Cp = vecComplex(p,N*Nb);
	
	//defining DDS and minusDS
	spMat DDS(2*N*Nb+1,2*N*Nb+1);
	vec minusDS(2*N*Nb+1);
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//beginning newton-raphson loop
	while ((sol_test.back()>closenessS || solM_test.back()>closenessSM || runs_count<min_runs))
		{
		runs_count ++;
		
		//defining the zero mode at the final time boundary and the time step before
		vec Chi0(Nb*N);
		Chi0 = Eigen::VectorXd::Zero(N*Nb);
		//#pragma omp parallel for
		for (unsigned int j=0; j<N; j++)
			{
			unsigned int pos = (j+1)*Nb-1;
			long int neighPos = neigh(pos,1,1,Nb,N);
			if(neighPos!=-1)
            	{
            	Chi0(pos) = p(2*neighPos)-p(2*neigh(pos,1,-1,Nb,N)); //final time slice
            	//Chi0(pos-1) = p(2*neigh(pos-1,1,1,Nb,N))-p(2*neigh(pos-1,1,-1,Nb,N)); //penultimate time slice
            	}
            else
            	{
            	Chi0(pos) = 0.0;
            	}
            }
        if (runs_count==1) //printing Chi0
        	{
        	//printVectorB("data/" + timeNumber + "Chi0.dat",Chi0);
        	}

		// allocating memory for DS, DDS
		minusDS = Eigen::VectorXd::Zero(2*N*Nb+1); //initializing to zero
		DDS.setZero(); //just making sure
		Eigen::VectorXi DDS_to_reserve(2*N*Nb+1);//number of non-zero elements per column
		DDS_to_reserve = Eigen::VectorXi::Constant(2*N*Nb+1,11);
		DDS_to_reserve(0) = 3; //these need to be changed when boundary conditions need to be more compicated
		DDS_to_reserve(1) = 3;
		DDS_to_reserve(2*N*Nb-2) = 3;
		DDS_to_reserve(2*N*Nb-1) = 3;
		DDS_to_reserve(2*N*Nb) = N;
		DDS.reserve(DDS_to_reserve);
		
		//initializing to zero
		comp kineticS = 0.0;
		comp kineticT = 0.0;
		comp pot_0 = 0.0;
		comp pot_r = 0.0;
		erg = Eigen::VectorXcd::Constant(NT,-ergZero);
		
		//testing that the potential term is working for pot3
		if (pot[0]=='3' && false)
			{
			comp Vtrial = 0.0, Vcontrol = 0.0;
			for (unsigned int j=0; j<N; j++)
				{
				double r = 1.0e-16 + j*a;
				paramsV  = {r, 0.0};
				Vcontrol += pow(p(2*j*Nb),2.0)/2.0 - pow(p(2*j*Nb),4.0)/4.0/pow(r,2.0);
				Vtrial += V(p(2*j*Nb));
				}
			double potTest = pow(pow(real(Vcontrol-Vtrial),2.0) + pow(imag(Vcontrol-Vtrial),2.0),0.5);
			cout << "potTest = " << potTest << endl;
			}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//assigning values to minusDS and DDS and evaluating action
		//#pragma omp parallel for
		for (unsigned long int j = 0; j < N*Nb; j++)
			{		
			unsigned int t = intCoord(j,0,Nb); //coordinates
			unsigned int x = intCoord(j,1,Nb);
			long int neighPosX = neigh(j,1,1,Nb,N);
			if (pot[0]=='3')
				{
				paramsV  = {1.0e-16+x*a, 0.0};
				}

			
			if (absolute(Chi0(j))>1.0e-16)  //zero mode lagrange constraint
				{
				DDS.insert(2*j,2*N*Nb) = a*Chi0(j); 
				DDS.insert(2*N*Nb,2*j) = a*Chi0(j);
		    	minusDS(2*j) += -a*Chi0(j)*p(2*N*Nb);
		    	minusDS(2*N*Nb) += -a*Chi0(j)*p(2*j);
		    	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//boundaries
			if (pot[0]=='3' && x==(N-1))
				{
				DDS.insert(2*j,2*j) = 1.0; //p=0 at r=R
				DDS.insert(2*j+1,2*j+1) = 1.0;
				}
			else if (pot[0]=='3' && x==0)
				{
				DDS.insert(2*j,2*j) = -1.0/a; //dp/dx=1 at r=0
				DDS.insert(2*j,2*(j+1)) = 1.0/a;
				DDS.insert(2*j+1,2*j+1) = -1.0/a;
				DDS.insert(2*j+1,2*(j+1)+1) = 1.0/a;
				}
			else if (t==(Nb-1))
				{
				comp Dt = -b*i/2.0;
				erg(t+Na) += pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				kineticS += Dt*pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0; //n.b. no contribution from time derivative term at the final time boundary	
				pot_0 += Dt*a*V(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				
				
				DDS.insert(2*j,2*j) = 1.0/b; //zero time derivative
				DDS.insert(2*j,2*(j-1)) = -1.0/b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
			else if (t==0)
				{
				comp dt = -b*i;
				comp Dt = -b*i/2.0;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				kineticS += Dt*pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0;
				pot_0 += Dt*a*V(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				erg(t+Na) += a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
				if (inP.compare("b")==0)
					{
					DDS.insert(2*j,2*j) = 1.0; //zero change (initial input satisfies b.c.s)
					DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
					}
				else if (inP.compare("p")==0)
					{
					DDS.insert(2*j,2*j) = -1.0/b; //zero time derivative
					DDS.insert(2*j,2*(j+1)) = 1.0/b;
					DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
					}
				}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//bulk
			else
				{
				comp dt = -b*i;
				comp Dt = -b*i;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				kineticS += Dt*pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0;
				pot_0 += Dt*a*V(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				erg(t+Na) += a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
                for (unsigned int k=0; k<2*2; k++)
                	{
                    int sign = pow(-1,k);
                    int direc = (int)(k/2.0);
                    long int neighb = neigh(j,direc,sign,Nb,N);
                    if (direc == 0)
                    	{
                        minusDS(2*j) += real(a*Cp(j+sign)/dt);
                        minusDS(2*j+1) += imag(a*Cp(j+sign)/dt);
                        DDS.insert(2*j,2*(j+sign)) = -real(a/dt);
                        DDS.insert(2*j,2*(j+sign)+1) = imag(a/dt);
                        DDS.insert(2*j+1,2*(j+sign)) = -imag(a/dt);
                        DDS.insert(2*j+1,2*(j+sign)+1) = -real(a/dt);
                        }
                    else if (neighb !=-1)
                    	{
                        minusDS(2*j) += - real(Dt*Cp(neighb)/a);
                        minusDS(2*j+1) += - imag(Dt*Cp(neighb)/a);
                        DDS.insert(2*j,2*neighb) = real(Dt/a);
                        DDS.insert(2*j,2*neighb+1) = -imag(Dt/a);
                        DDS.insert(2*j+1,2*neighb) = imag(Dt/a);
                        DDS.insert(2*j+1,2*neighb+1) = real(Dt/a);
                        }
                    }
                comp temp0 = 2.0*a/dt;
	            comp temp1 = a*Dt*(2.0*Cp(j)/pow(a,2.0) + dV(Cp(j)) + dVr(Cp(j)));
	            comp temp2 = a*Dt*(2.0/pow(a,2.0) + ddV(Cp(j)) + ddVr(Cp(j)));
	                
	            minusDS(2*j) += real(temp1 - temp0*Cp(j));
	            minusDS(2*j+1) += imag(temp1 - temp0*Cp(j));
	            DDS.insert(2*j,2*j) = real(-temp2 + temp0);
	            DDS.insert(2*j,2*j+1) = imag(temp2 - temp0);
	            DDS.insert(2*j+1,2*j) = imag(-temp2 + temp0);
	            DDS.insert(2*j+1,2*j+1) = real(-temp2 + temp0);
	            }
            }
        action = kineticT - kineticS - pot_0 - pot_r;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			string prefix = "./data/" + timeNumber;
			string suffix = inP+"_" + numberToString<unsigned int>(loop)+"_" + numberToString<unsigned int>(runs_count)+".dat";
			if ((print_choice.compare("v")==0 || print_choice.compare("e")==0))
				{
				string minusDSfile = prefix + "minusDSE"+suffix;
				printVectorB(minusDSfile,minusDS);
				}
			if ((print_choice.compare("p")==0 || print_choice.compare("e")==0))
				{
				string piEarlyFile = prefix + "piE"+suffix;
				printVectorB(piEarlyFile,p);
				}
			if ((print_choice.compare("m")==0 || print_choice.compare("e")==0))
				{
				string DDSfile = prefix + "DDSE"+suffix;
				printSpmat(DDSfile,DDS);
				}
			}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

	//solving for delta in DDS*delta=minusDS, where p' = p + delta		
		vec delta(2*N*Nb+1);
		delta = Eigen::VectorXd::Zero(2*N*Nb+1);
		DDS.makeCompressed();
		Eigen::SparseLU<spMat> solver;
		
		solver.analyzePattern(DDS);
		if(solver.info()!=Eigen::Success)
			{
			cerr << "DDS pattern analysis failed, solver.info() = "<< solver.info() << endl;
			return 0;
			}		
		solver.factorize(DDS);
		if(solver.info()!=Eigen::Success) 
			{
			cerr << "Factorization failed, solver.info() = "<< solver.info() << endl;
			return 0;
			}
		delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side
		if(solver.info()!=Eigen::Success)
			{
			cerr << "Solving failed, solver.info() = "<< solver.info() << endl;
			cerr << "log(abs(det(DDS))) = " << solver.logAbsDeterminant() << endl;
			cerr << "sign(det(DDS)) = " << solver.signDeterminant() << endl;
			return 0;
			}
		
		//independent check on whether calculation worked
		vec diff(2*N*Nb+1);
		diff = DDS*delta-minusDS;
		double maxDiff = diff.maxCoeff();
		maxDiff = absolute(maxDiff);
		calc_test.push_back(maxDiff);
		if (calc_test.back()>closenessC)
			{
			cout << "Calculation failed" << endl;
			cout << "calc_test = " << calc_test.back() << endl;
			return 0;
			}

		//assigning values to phi
		p += delta;
		
		//passing changes on to complex vector
		Cp = vecComplex(p,N*Nb);
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//convergence issues
		
		//checking pot_r is much smaller than the other potential terms
		reg_test.push_back(abs(pot_r/pot_0));
		if (reg_test.back()>closenessR)
			{
			cout << "regularisation term is too large, regTest = " << reg_test.back() << endl;
			}
		
		//evaluating norms
		double normDS = minusDS.dot(minusDS);
		normDS = pow(normDS,0.5);
		double maxDS = minusDS.maxCoeff();
		double minDS = minusDS.minCoeff();
		if (-minDS>maxDS)
			{
			maxDS = -minDS;
			}
		double normP = p.dot(p);
		normP = pow(normP,0.5);
		double normDelta = delta.dot(delta);
		normDelta = pow(normDelta,0.5);
		
		//assigning test values
		//quantities used to stop newton-raphson loop
		action_test.push_back(abs(action - action_last)/abs(action_last));
		action_last = action;
		sol_test.push_back(normDS/normP);
		solM_test.push_back(maxDS);
		delta_test.push_back(normDelta/normP);
			
		//printing tests to see convergence
		if (runs_count==1)
			{
			printf("%16s%16s%16s%16s%16s%16s\n","loop","runsCount","actionTest","solTest","solMTest","deltaTest");
			}
		printf("%16i%16i%16g%16g%16g%16g\n",loop,runs_count,action_test.back(),sol_test.back(),solM_test.back(),delta_test.back());
		
		} //closing "runs" while loop
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//propagating solution along minkowskian time paths
	
	//propagating solution back in minkowskian time
    //A1. initialize mp==mphi using last point of ephi and zeros- use complex phi
    cVec ap(N*(Na+1)); //phi on section "a"
    ap = Eigen::VectorXcd::Zero(N*(Na+1));
    //#pragma omp parallel for
    for (unsigned int j=0; j<N; j++)
    	{
        ap(j*(Na+1)) = Cp(j*Nb);
    	}

    //A2. initialize vel - defined at half steps, first step being at t=-1/2,
    //vel(t+1/2) := (p(t+1)-p(t))/dt
    cVec velA (N*(Na+1));
    velA = Eigen::VectorXcd::Zero(N*(Na+1));
    double dtau = -b;
    double Dt0 = dtau; //b/2*(-1+1i); - this is surely wrong!!
    //#pragma omp parallel for
    //for (unsigned int j=0; j<N; j++)
    	//{
        //velA(j*(Na+1)) = 0; //due to boundary condition
    	//}

    
    //A3. initialize acc using phi and expression from equation of motion and zeros-complex
    cVec accA(N*(Na+1));
    accA = Eigen::VectorXcd::Zero(N*(Na+1));
    //#pragma omp parallel for
    for (unsigned int j=0; j<N; j++)
    	{
    	if (pot[0]=='3')
			{
			paramsV  = {1.0e-16+j*a, A};
			}
		unsigned int l = j*(Na+1);
		if (pot[0]=='3' && j==(N-1))
			{
			accA(l) = 0.0;
			}
		else if (pot[0]=='3' && j==0)
			{
			accA(l) = ((Dt0/pow(a,2.0))*(2.0*ap(neigh(l,1,1,Na+1,N))-2.0*ap(l)) \
            	-Dt0*(dV(ap(l))+dVr(ap(l))))/dtau;
			}
    	else
    		{
        	accA(l) = ((Dt0/pow(a,2.0))*(ap(neigh(l,1,1,Na+1,N))+ap(neigh(l,1,-1,Na+1,N))-2.0*ap(l)) \
            	-Dt0*(dV(ap(l))+dVr(ap(l))))/dtau;
            }
    	}
    	
    //A4.5 starting the energy and that off
	vec linErgA(Na); linErgA = Eigen::VectorXd::Zero(Na);
	vec linNumA(Na); linNumA = Eigen::VectorXd::Zero(Na);

    //A7. run loop
    for (unsigned int t=1; t<(Na+1); t++)
    	{
    	//#pragma omp parallel for
        for (unsigned int x=0; x<N; x++)
        	{
            unsigned int m = t+x*(Na+1);
            velA(m) = velA(m-1) + dtau*accA(m-1);
            ap(m) = ap(m-1) + dtau*velA(m);
        	}
        //#pragma omp parallel for	
        for (unsigned int x=0; x<N; x++)
        	{
        	if (pot[0]=='3')
				{
				paramsV  = {1.0e-16+x*a, A};
				}
            unsigned int m = t+x*(Na+1);
            if (pot[0]=='3' && x==(N-1))
				{
				accA(m) = 0.0;
				}
			else if (pot[0]=='3' && x==0)
				{
				accA(m) = (1.0/pow(a,2.0))*(2.0*ap(neigh(m,1,1,Na+1,N))-2.0*ap(m)) \
            		-dV(ap(m)) - dVr(ap(m));
            	erg (Na-t) += a*pow(ap(m-1)-ap(m),2.0)/pow((-dtau),2.0)/2.0 + pow(ap(neigh(m,1,1,Na+1,N))-ap(m),2.0)/a/2.0 \
            		+ a*V(ap(m)) + a*Vr(ap(m));
				}
			else
				{
		    	accA(m) = (1.0/pow(a,2.0))*(ap(neigh(m,1,1,Na+1,N))+ap(neigh(m,1,-1,Na+1,N))-2.0*ap(m)) \
            		-dV(ap(m)) - dVr(ap(m));
            	erg (Na-t) += a*pow(ap(m-1)-ap(m),2.0)/pow((-dtau),2.0)/2.0 + pow(ap(neigh(m,1,1,Na+1,N))-ap(m),2.0)/a/2.0 \
            		+ a*V(ap(m)) + a*Vr(ap(m));
		        }
            for (unsigned int y=0; y<N; y++)
            	{
            	unsigned int n = t + y*(Na+1);
		        if (absolute(theta)<1.0e-16)
					{
		        	linErgA(Na-t) += Eomega(x,y)*(real(ap(m))-minima[0])*(real(ap(n))-minima[0]) + Eomega(x,y)*imag(ap(m))*imag(ap(n));
					linNumA (Na-t) += omega(x,y)*(real(ap(m))-minima[0])*(real(ap(n))-minima[0]) + omega(x,y)*imag(ap(m))*imag(ap(n));
		        	}
				else
					{
					linErgA(Na-t) += 2.0*Gamma*omega(x,y)*(real(ap(m))-minima[0])*(real(ap(n))-minima[0])/pow(1.0+Gamma,2.0) + 2.0*Gamma*omega(x,y)*imag(ap(m))*imag(ap(n))/pow(1.0-Gamma,2.0);
					linNumA(Na-t) += 2.0*Gamma*Eomega(x,y)*(real(ap(m))-minima[0])*(real(ap(n))-minima[0])/pow(1.0+Gamma,2.0) + 2.0*Gamma*Eomega(x,y)*imag(ap(m))*imag(ap(n))/pow(1.0-Gamma,2.0);
					}
				}
        	}
    	}
		
	//E = 0;
	//unsigned int linearInt = (unsigned int)(Na/6);
	//#pragma omp parallel for
	//for (unsigned int j=0; j<linearInt; j++)
	//	{
	//	E += real(erg(j));
	//	}
	//E /= linearInt;
	E = linErgA(0);
	W = - E*2.0*Tb + 2.0*imag(action);

    //now propagating forwards along c
    //C2. initialize mp==mphi using last point of ephi and zeros- use complex phi
    cVec ccp(N*(Nc+1)); //phi on section "c"
    ccp = Eigen::VectorXcd::Zero(N*(Nc+1));
    //#pragma omp parallel for
    for (unsigned int j=0; j<N; j++)
    	{
        ccp(j*(Nc+1)) = Cp(j*Nb+Nb-1);
    	}

    //C3. initialize vel - defined at half steps, first step being at t=-1/2,
    //vel(t+1/2) := (p(t+1)-p(t))/dt
    cVec velC (N*(Nc+1));
    velC = Eigen::VectorXcd::Zero(N*(Nc+1));
    dtau = b;
    Dt0 = dtau; //b/2*(1-1i); - this is surely wrong!!
    //#pragma omp parallel for
    //for (unsigned int j=0; j<N; j++)
    	//{
        //velC(j*(Nc+1)) = 0; //due to boundary condition
    	//}

    //C4. initialize acc using phi and expression from equation of motion and zeros-complex
    cVec accC(N*(Nc+1));
    accC = Eigen::VectorXcd::Zero(N*(Nc+1));
    //#pragma omp parallel for
    for (unsigned int j=0; j<N; j++)
    	{
    	unsigned int l = j*(Nc+1);
    	if (pot[0]=='3')
			{
			paramsV  = {1.0e-16+j*a, A};
			}
		if (pot[0]=='3' && j==(N-1))
			{
			accC(l) = 0.0;
			}
		else if (pot[0]=='3' && j==0)
			{
			accC(l) = ((Dt0/pow(a,2.0))*(2.0*ccp(neigh(l,1,1,Nc+1,N))-2.0*ccp(l)) \
            		-Dt0*(dV(ccp(l))+dVr(ccp(l))))/dtau;
			}
    	else
    		{
        	accC(l) = ((Dt0/pow(a,2.0))*(ccp(neigh(l,1,1,Nc+1,N))+ccp(neigh(l,1,-1,Nc+1,N))-2.0*ccp(l)) \
            		-Dt0*(dV(ccp(l))+dVr(ccp(l))))/dtau;
            }
    	}

    //C7. run loop
    for (unsigned int t=1; t<(Nc+1); t++)
		{
		//#pragma omp parallel for
		for (unsigned int x=0; x<N; x++)
			{
		    unsigned int l = t+x*(Nc+1);
		    velC(l) = velC(l-1) + dtau*accC(l-1);
		    ccp(l) = ccp(l-1) + dtau*velC(l);
			}
		//#pragma omp parallel for
		for (unsigned int x=0; x<N; x++)
			{
			if (pot[0]=='3')
				{
				paramsV  = {1.0e-16+x*a, A};
				}
		    unsigned int l = t+x*(Nc+1);
		    if (pot[0]=='3' && x==(N-1))
				{
				accC(l) = 0.0;
				}
			else if (pot[0]=='3' && x==0)
				{
				accC(l) = (1.0/pow(a,2.0))*(2.0*ccp(neigh(l,1,1,Nc+1,N))-2.0*ccp(l)) \
		    		-dV(ccp(l));
				}
			else
				{
		    	accC(l) = (1.0/pow(a,2.0))*(ccp(neigh(l,1,1,Nc+1,N))+ccp(neigh(l,1,-1,Nc+1,N))-2.0*ccp(l)) \
		    		-dV(ccp(l));
		        }
		    if ((t>1) && x!=(N-1))
		    	{
				erg (Na+Nb-2+t) += a*pow(ccp(l)-ccp(l-1),2.0)/pow(dtau,2.0)/2.0\
				 	+ pow(ccp(neigh(l-1,1,1,Nc+1,N))-ccp(l-1),2.0)/a/2.0\
	    			+ a*V(ccp(l-1)) + a*Vr(ccp(l-1));
		    	}
			}
		}
		
	//checking energy conserved
	double ergChange = 0.0;
	double relErgChange = 0.0;
	if (absolute(real(erg(0)))>1.0e-16)
		{
		ergChange = absolute(real(erg(0))-real(erg(NT-2)));
		relErgChange = absolute((real(erg(0))-real(erg(NT-2)))/real(erg(0)));
		}
	erg_test.push_back(ergChange);
	if (erg_test.back()>closenessE)
		{
		cout << "energy change = " << ergChange << endl;
		cout << "relative energy change = " << relErgChange << endl;
		}
    

    //12. combine phi with ap and cp and save combination to file
    cVec tCp(NT*N);
    //#pragma omp parallel for
    for (unsigned int j=0; j<NT*N; j++)
    	{
        unsigned int t = intCoord(j,0,NT);
        unsigned int x = intCoord(j,1,NT);
        if (t<Na)
        	{
            t = Na-t;
            tCp(j) = ap(t+x*(Na+1));
            }
        else if (t<(Na+Nb))
        	{
            t = t - Na;
            tCp(j) = Cp(t+x*Nb);
            }
        else
        	{
            t = t - Na - Nb + 1;
            tCp(j) = ccp(t+x*(Nc+1));
        	}
    	}
    	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	    	//misc end of program tasks - mostly printing
    
    //making real vec from complex one
    vec tp(2*N*NT);
    tp = vecReal(tCp,NT*N);
    tp.conservativeResize(2*N*NT+1);
    tp(2*N*NT) = p(2*N*Nb);
    
    //stopping clock
	time = clock() - time;
	double realtime = time/1000000.0;
	
	//printing to terminal
	printf("\n");
	printf("%8s%8s%8s%8s%8s%8s%8s%16s%16s\n","runs","time","N","NT","L","Tb","dE","E","W");
	printf("%8i%8g%8i%8i%8g%8g%8g%16g%16g\n",runs_count,realtime,N,NT,L,Tb,dE,E,W);
	printf("\n");
	 printf("%60s\n","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

	//printing action value
	FILE * actionfile;
	actionfile = fopen("./data/action.dat","a");
	fprintf(actionfile,"%16s%8i%8i%8g%8g%8g%14g%14g%14g%14g\n",timeNumber.c_str(),N,NT,L,Tb,dE,E\
	,W,sol_test.back(),erg_test.back());
	fclose(actionfile);
	
	string prefix = "./data/" + timeNumber;
	string suffix = "_" + numberToString<unsigned int>(loop) + ".dat";
	
	//copying a version of inputs with timeNumber
	string runInputs = prefix + "inputsPi" + "_" + numberToString<unsigned int>(loop);
	if (loop_choice[0]=='N')
		{
		changeInputs(runInputs,loop_choice, numberToString<unsigned int>(intLoopParameter));
		}
	else if (loop_choice.compare("n")!=0)
		{
		changeInputs(runInputs,loop_choice, numberToString<double>(doubleLoopParameter));
		}
	else
		{
		copyFile("inputs",runInputs);
		}
	printf("%12s%30s\n","output: ",runInputs.c_str());

	//printing output phi on Euclidean time part
	string pifile = prefix + "pi"+inP+suffix;
	printVectorB(pifile,p);
	printf("%12s%30s\n"," ",pifile.c_str());
	
	//printing output phi on whole time contour
	string tpifile = prefix + "tpi"+inP+suffix;
	printVector(tpifile,tp);
	printf("%12s%30s\n"," ",tpifile.c_str());
	//gp(tpifile,"repi.gp");
	
	//printing output minusDS				
	string minusDSfile = "./data/" + timeNumber + "minusDS"+suffix;
	printVectorB(minusDSfile,minusDS);
	printf("%12s%30s\n"," ",minusDSfile.c_str());
				
	//printing output DDS
	string DDSfile = prefix + "DDS"+inP+suffix;
	printSpmat(DDSfile,DDS);
	printf("%12s%30s\n"," ",DDSfile.c_str());
	
	//printing linErgVec
	string linErgFile = "./data/" + timeNumber + "linErg"+inP+suffix;
	simplePrintVector(linErgFile,linErgA);
	printf("%12s%30s\n"," ",linErgFile.c_str());
//	gpSimple(linErgFile);
	
	//printing erg
	string ergFile = prefix + "erg"+inP+suffix;
	simplePrintCVector(ergFile,erg);
	printf("%12s%30s\n"," ",ergFile.c_str());
//	gpSimple(ergFile);
	
	//printing error, and eigenvalue to file
	ifstream eigenvalueIn;
	string eigenvaluefile = "data/eigValue.dat";
	eigenvalueIn.open(eigenvaluefile.c_str(),ios::in);
	string lastEigLine = getLastLine(eigenvalueIn);
	eigenvalueIn.close();
	ofstream eigenvalueOut;
	string eigValueFile = prefix + "eigValue.dat";
	eigenvalueOut.open(eigValueFile.c_str(),ios::out);
	eigenvalueOut << lastEigLine;
	eigenvalueOut.close();
	printf("%12s%30s\n"," ",eigenvaluefile.c_str());

	//printing eigenvector to file
	string eigenvectorFile = prefix + "eigVec.dat";
	copyFile("data/eigVec.dat",eigenvectorFile);
	printf("%12s%30s\n"," ",eigenvectorFile.c_str());

} //closing parameter loop

return 0;
}
