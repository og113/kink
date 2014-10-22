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
string timeNumber = currentDateTime();

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
			if (absolute(theta)>2.0e-16)
				{
				cout << "theta != 0" << endl;
				cout << "theta = " << theta << endl;
				}
			}
		else if (lineNumber==1)
			{
			ss >> aq.inputChoice >> aq.inputFile >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
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
	
//determining number of runs
closenessA = 1.0;
closenessS = 1.0e-5;
closenessSM = 1.0e-4;
closenessD = 1.0;
closenessC = 1.0e-12;
closenessE = 1.0e-2;
closenessR = 1.0e-4;

string loop_choice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
string print_choice = aq.printChoice;
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//calculated quantities derived from the inputs, and loading eigVec and eigVal	

//derived quantities
NT = Na + Nb + Nc;
Gamma = exp(-theta);
epsilon = dE;
R = 2.0/3.0/epsilon;
alpha *= R;
L = LoR*R;
vec negVec(2*N*Nb+1);
unsigned int negP;
double negc;
double negcheck;
double negerror; //should be <<1
if (inP.compare("p") == 0 || inP.compare("f") == 0)
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
			system("./negEig");
			system(timeNumber.c_str()); //won't work
			cout << "negEig run" << endl;
			}
		negVec = loadVector("./data/eigVec.dat",Nb,1);
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
else if (inP.compare("b") == 0)
	{
	Tb = 1.5*R;
	}
a = L/(N-1.0);
b = Tb/(Nb-1.0);
Ta = b*Na;
Tc = b*Nc;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//begin loop over varying parameter
//(NB. may only be one loop)
//and defining some quantities to be used later
for (unsigned int loop=0; loop<aq.totalLoops; loop++)
	{
	//giving values of varying parameters
	int intLoopParameter;
	double doubleLoopParameter;
	if (loop_choice.compare("N") == 0)
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
	
	//potential functions
	if (pot[0]=='1')
		{
		V_params = &V1_params;
		V = &V1;
		V0 = &V10;
		Ve = &V1e;
		dV_params = &dV1_params;
		dV = &dV1;
		ddV_params = &ddV1_params;
		ddV = &ddV1;
		}
	else if (pot[0]=='2')
		{
		V_params = &V2_params;
		V = &V2;
		V0 = &V20;
		Ve = &V2e;
		dV_params = &dV2_params;
		dV = &dV2;
		ddV_params = &ddV2_params;
		ddV = &ddV2;
		}
	
	//defining some important scalar quantities
	double S1 = 2.0/3.0; //mass of kink multiplied by lambda
	double twaction = -pi*epsilon*pow(R,2)/2.0 + pi*R*S1;
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
	
	//finding roots of dV
	struct f_gsl_params params = { epsilon, A};
	gsl_function_fdf FDF;
	FDF.f = f_gsl;
	FDF.df = df_gsl;
	FDF.fdf = fdf_gsl;
	FDF.params = &params;
	
	root = minimaFn(&FDF, -3.0, 3.0, 20);
	sort(root.begin(),root.end());
	comp ergZero = N*a*V(root[0]);
	mass2 = real(ddV(root[0]));
	
	//defining lambda functions for regularization
	auto Vr = [&] (const comp & phi)
		{
		return -i*reg*VrFn(phi,root[0],root[2]);
		};
	auto dVr = [&] (const comp & phi)
		{
		return -i*reg*dVrFn(phi,root[0],root[2]);
		};
	auto ddVr = [&] (const comp & phi)
		{
		return -i*reg*ddVrFn(phi,root[0],root[2]);
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
		
	for (unsigned int j=0; j<N; j++)
		{
		for (unsigned int k=0; k<N; k++)
			{
			for (unsigned int l=0; l<N; l++)
				{
				omega(j,k) += a*pow(eigenValues(l),0.5)*eigenVectors(j,l)*eigenVectors(k,l);
				Eomega(j,k) += a*eigenValues(l)*eigenVectors(j,l)*eigenVectors(k,l);
				}
			}
		}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//assigning input phi
	if (inP.compare("f")!=0 && loop==0)
		{
		if (R<alpha)
			{
			cout << "R is too small. Not possible to give thinwall input. It should be more that " << alpha;
			}
		for (unsigned int j=0; j<N*Nb; j++)
			{
			comp t = coordB(j,0);
			comp x = coordB(j,1);
			p(2*j+1) = 0.0; //imaginary parts set to zero
			if (inP.compare("b")==0 || (inP.compare("p")==0 && Tb>R))
				{
				double rho = real(sqrt(-pow(t,2.0) + pow(x,2.0))); //should be real even without real()
				if ((rho-R)<-alpha)
					{
					p(2*j) = root[2];
					}
				else if ((rho-R)>alpha)
					{
					p(2*j) = root[0];
					}
				else
					{
					p(2*j) = (root[2]+root[0])/2.0 + (root[0]-root[2])*tanh((rho-R)/2.0)/2.0;
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
					p(2*j) = root[2];
					}
				else if ((rho1-R)>alpha || (rho2-R)>alpha)
					{
					p(2*j) = root[0];
					}
				else if (real(x)>0) //note that the coord should be real
					{
					p(2*j) = (root[2]+root[0])/2.0 + (root[0]-root[2])*tanh((rho1-R)/2.0)/2.0;
					}
				else if (real(x)<0)
					{
					p(2*j) = (root[2]+root[0])/2.0 + (root[0]-root[2])*tanh((rho2-R)/2.0)/2.0;
					}
				else
					{
					p(2*j) = root[2]; //i.e. if coordB(j,1) == 0
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
			loadfile = "./data/" + aq.inputFile + ".dat";
			inP = "p";
			}
		else
			{
			loadfile = "./data/" + timeNumber + "pi"+inP+to_string(loop-1)+".dat";
			}
		
		p = loadVector(loadfile,Nb,1);
		}
	
	//fixing input periodic instanton to have zero time derivative at time boundaries
    if (true)
    	{
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
		
	//very early vector print
	string earlyPrintFile = "data/" + timeNumber + "piE"+inP+ to_string(loop) + "0.dat";
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
		for (unsigned int j=0; j<N; j++)
			{
			unsigned int pos = (j+1)*Nb-1;
            Chi0(pos) = p(2*neigh(pos,1,1,Nb))-p(2*neigh(pos,1,-1,Nb)); //final time slice
            //Chi0(pos-1) = p(2*neigh(pos-1,1,1,Nb))-p(2*neigh(pos-1,1,-1,Nb)); //penultimate time slice
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
		comp pot_l = 0.0;
		comp pot_e = 0.0;
		comp pot_r = 0.0;
		erg = Eigen::VectorXcd::Constant(NT,-ergZero);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//assigning values to minusDS and DDS and evaluating action
		for (unsigned long int j = 0; j < N*Nb; j++)
			{		
			unsigned int t = intCoord(j,0,Nb); //coordinates
			
			if (absolute(Chi0(j))>2.0e-16) //zero mode lagrange constraint
				{
				DDS.insert(2*j,2*N*Nb) = a*Chi0(j); 
				DDS.insert(2*N*Nb,2*j) = a*Chi0(j);
		    	minusDS(2*j) += -a*Chi0(j)*p(2*N*Nb);
		    	minusDS(2*N*Nb) += -a*Chi0(j)*p(2*j);
		    	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//boundaries			
			if (t==(Nb-1))
				{
				comp Dt = -b*i/2.0;
				kineticS += Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0; //n.b. no contribution from time derivative term at the final time boundary
				pot_l += Dt*a*V0(Cp(j));
				pot_e += Dt*a*Ve(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				erg(t+Na) += pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
				DDS.insert(2*j,2*j) = 1.0/b; //zero time derivative
				DDS.insert(2*j,2*(j-1)) = -1.0/b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
			else if (t==0)
				{
				comp dt = -b*i;
				comp Dt = -b*i/2.0;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				kineticS += Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += Dt*a*V0(Cp(j));
				pot_e += Dt*a*Ve(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				erg(t+Na) += a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
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
				kineticS += Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += Dt*a*V0(Cp(j));
				pot_e += Dt*a*Ve(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				erg(t+Na) += a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
                for (unsigned int k=0; k<2*2; k++)
                	{
                    int sign = pow(-1,k);
                    int direc = (int)(k/2.0);
                    if (direc == 0)
                    	{
                        minusDS(2*j) += real(a*Cp(j+sign)/dt);
                        minusDS(2*j+1) += imag(a*Cp(j+sign)/dt);
                        DDS.insert(2*j,2*(j+sign)) = -real(a/dt);
                        DDS.insert(2*j,2*(j+sign)+1) = imag(a/dt);
                        DDS.insert(2*j+1,2*(j+sign)) = -imag(a/dt);
                        DDS.insert(2*j+1,2*(j+sign)+1) = -real(a/dt);
                        }
                    else
                    	{
                        unsigned int neighb = neigh(j,direc,sign,Nb);
                        
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
        action = kineticT - kineticS - pot_l - pot_e - pot_r;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			string prefix = "./data/" + timeNumber;
			string suffix = to_string(loop)+"_" + to_string(runs_count)+".dat";
			if ((print_choice.compare("a")==0 || print_choice.compare("e")==0) && 1==0) //have stopped this one as it's annoying
				{
				printAction(kineticT-kineticS,pot_l,pot_e);
				}
			if ((print_choice.compare("v")==0 || print_choice.compare("e")==0) && 1==0)
				{
				string minusDSfile = prefix + "minusDSE"+inP+suffix;
				printVectorB(minusDSfile,minusDS);
				}
			if ((print_choice.compare("p")==0 || print_choice.compare("e")==0) && delta_test.back()>0.2 && 1==0)
				{
				string piEarlyFile = prefix + "piE"+inP+to_string(loop)+suffix;
				printVectorB(piEarlyFile,p);
				}
			if ((print_choice.compare("m")==0 || print_choice.compare("e")==0) && 1==0)
				{
				string DDSfile = prefix + "DDSE"+inP+to_string(loop)+suffix;
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
			cout << "DDS pattern analysis failed, solver.info() = "<< solver.info() << endl;
			return 0;
			}		
		solver.factorize(DDS);
		if(solver.info()!=Eigen::Success) 
			{
			cout << "Factorization failed, solver.info() = "<< solver.info() << endl;
			return 0;
			}
		delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side
		if(solver.info()!=Eigen::Success)
			{
			cout << "Solving failed, solver.info() = "<< solver.info() << endl;
			cout << "log(abs(det(DDS))) = " << solver.logAbsDeterminant() << endl;
			cout << "sign(det(DDS)) = " << solver.signDeterminant() << endl;
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
		reg_test.push_back(abs(pot_r/pot_l));
		if (reg_test.back()<abs(pot_r/pot_e))
			{
			reg_test.back() = abs(pot_r/pot_e);
			}

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
    for (unsigned int j=0; j<N; j++)
    	{
        velA(j*(Na+1)) = 0; //due to boundary condition
    	}

    
    //A3. initialize acc using phi and expression from equation of motion and zeros-complex
    cVec accA(N*(Na+1));
    accA = Eigen::VectorXcd::Zero(N*(Na+1));
    for (unsigned int j=0; j<N; j++)
    	{
    	unsigned int l = j*(Na+1);
        accA(l) = ((Dt0/pow(a,2.0))*(ap(neigh(l,1,1,Na+1))+ap(neigh(l,1,-1,Na+1))-2.0*ap(l)) \
            -Dt0*(dV(ap(l))+dVr(ap(l))))/dtau;
    	}
    	
    //A4.5 starting the energy and that off
	vec linErgA(Na); linErgA = Eigen::VectorXd::Zero(Na);
	vec linNumA(Na); linNumA = Eigen::VectorXd::Zero(Na);

    //A7. run loop
    for (unsigned int t=1; t<(Na+1); t++)
    	{
        for (unsigned int x=0; x<N; x++)
        	{
            unsigned int m = t+x*(Na+1);
            velA(m) = velA(m-1) + dtau*accA(m-1);
            ap(m) = ap(m-1) + dtau*velA(m);
        	}
        for (unsigned int x=0; x<N; x++)
        	{
            unsigned int m = t+x*(Na+1);
            accA(m) = (1.0/pow(a,2.0))*(ap(neigh(m,1,1,Na+1))+ap(neigh(m,1,-1,Na+1))-2.0*ap(m)) \
            -dV(ap(m)) - dVr(ap(m));
            erg (Na-t) += a*pow(ap(m-1)-ap(m),2.0)/pow((-dtau),2.0)/2.0 + pow(ap(neigh(m,1,1,Na+1))-ap(m),2.0)/a/2.0 \
            + a*V(ap(m)) + a*Vr(ap(m));
            for (unsigned int y=0; y<N; y++)
            	{
            	unsigned int n = t + y*(Na+1);
		        if (absolute(theta)<2.0e-16)
					{
		        	linErgA(Na-t) += Eomega(x,y)*(real(ap(m))-root[0])*(real(ap(n))-root[0]) + Eomega(x,y)*imag(ap(m))*imag(ap(n));
					linNumA (Na-t) += omega(x,y)*(real(ap(m))-root[0])*(real(ap(n))-root[0]) + omega(x,y)*imag(ap(m))*imag(ap(n));
		        	}
				else
					{
					linErgA(Na-t) += 2.0*Gamma*omega(x,y)*(real(ap(m))-root[0])*(real(ap(n))-root[0])/pow(1.0+Gamma,2.0) + 2.0*Gamma*omega(x,y)*imag(ap(m))*imag(ap(n))/pow(1.0-Gamma,2.0);
					linNumA(Na-t) += 2.0*Gamma*Eomega(x,y)*(real(ap(m))-root[0])*(real(ap(n))-root[0])/pow(1.0+Gamma,2.0) + 2.0*Gamma*Eomega(x,y)*imag(ap(m))*imag(ap(n))/pow(1.0-Gamma,2.0);
					}
				}
        	}
    	}
		
	E = 0;
	unsigned int linearInt = (unsigned int)(Na/6);
	for (unsigned int j=0; j<linearInt; j++)
		{
		E += real(linErgA(j));
		}
	E /= linearInt;
	W = - E*2.0*Tb + 2.0*imag(action);

    //now propagating forwards along c
    //C2. initialize mp==mphi using last point of ephi and zeros- use complex phi
    cVec ccp(N*(Nc+1)); //phi on section "c"
    ccp = Eigen::VectorXcd::Zero(N*(Nc+1));
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
    for (unsigned int j=0; j<N; j++)
    	{
        velC(j*(Nc+1)) = 0; //due to boundary condition
    	}

    //C4. initialize acc using phi and expression from equation of motion and zeros-complex
    cVec accC(N*(Nc+1));
    accC = Eigen::VectorXcd::Zero(N*(Nc+1));
    for (unsigned int j=0; j<N; j++)
    	{
    	unsigned int l = j*(Nc+1);
        accC(l) = ((Dt0/pow(a,2.0))*(ccp(neigh(l,1,1,Nc+1))+ccp(neigh(l,1,-1,Nc+1))-2.0*ccp(l)) \
            -Dt0*(dV(ccp(l))+dVr(ccp(l))))/dtau;
    	}

    //C7. run loop
    for (unsigned int t=1; t<(Nc+1); t++)
		{
		for (unsigned int x=0; x<N; x++)
			{
		    unsigned int l = t+x*(Nc+1);
		    velC(l) = velC(l-1) + dtau*accC(l-1);
		    ccp(l) = ccp(l-1) + dtau*velC(l);
			}
		for (unsigned int x=0; x<N; x++)
			{
		    unsigned int l = t+x*(Nc+1);
		    accC(l) = (1.0/pow(a,2.0))*(ccp(neigh(l,1,1,Nc+1))+ccp(neigh(l,1,-1,Nc+1))-2.0*ccp(l)) \
		    -dV(ccp(l));
		    if (t>1)
		    	{
		    	erg (Na+Nb-2+t) += a*pow(ccp(l)-ccp(l-1),2.0)/pow(dtau,2.0)/2.0 + pow(ccp(neigh(l-1,1,1,Nc+1))-ccp(l-1),2.0)/a/2.0\
		    	+ a*V(ccp(l-1)) + a*Vr(ccp(l-1));
		    	}
			}
		}
		
	//checking energy conserved
	double ergChange = 0.0;
	double relErgChange = 0.0;
	if (absolute(real(erg(0)))>2.0e-16)
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
	string suffix = "_" + to_string(loop) + ".dat";
	
	//copying a version of inputs with timeNumber
	string runInputs = prefix + "inputsPi" + "_" + to_string(loop);
	if (loop_choice.compare("N") == 0)
		{
		changeInputs(runInputs,"N", to_string(intLoopParameter));
		}
	else if (loop_choice.compare("n")!=0)
		{
		changeInputs(runInputs,loop_choice, to_string(doubleLoopParameter));
		}
	else
		{
		copyFile("inputs",runInputs);
		}

	//printing output phi on Euclidean time part
	string pifile = prefix + "pi"+inP+suffix;
	printVectorB(pifile,p);
	
	//printing output phi on whole time contour
	string tpifile = prefix + "tpi"+inP+suffix;
	printVector(tpifile,tp);
	gp(tpifile,"repi.gp");
	
	//printing output minusDS				
	//string minusDSfile = "./data/" + timeNumber + "minusDS"+inP+to_string(loop)+".dat";
	//printVectorB(minusDSfile,minusDS);
				
	//printing output DDS
	string DDSfile = prefix + "DDS"+inP+suffix;
	printSpmat(DDSfile,DDS);
	
	//printing linErgVec
	string linErgFile = "./data/" + timeNumber + "linErg"+inP+suffix;
	simplePrintVector(linErgFile,linErgA);
//	gpSimple(linErgFile);
	
	//printing erg
	string ergFile = prefix + "erg"+inP+suffix;
	simplePrintCVector(ergFile,erg);
//	gpSimple(ergFile);
	
	//printing error, and eigenvalue to file
	FILE * eigenvalueFile;
	string eigValueFile = prefix + "eigValue.dat";
	eigenvalueFile = fopen(eigValueFile.c_str(),"w");
	fprintf(eigenvalueFile,"%16i%16g%16g%16g%16g\n",negP,negc,negcheck,negerror,negVal);
	fclose(eigenvalueFile);

	//printing eigenvector to file
	string eigenvectorFile = prefix + "eigVec.dat";
	printVectorB(eigenvectorFile,negVec);

} //closing parameter loop

return 0;
}
