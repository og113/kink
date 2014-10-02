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
#include <gsl/gsl_poly.h>
#include "gnuplot_i.hpp"
#include "pf.h"

using namespace std;

int main()
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//getting variables and user inputs from inputs

ifstream fin;
fin.open("inputs", ios::in);
string line;
int firstLine = 0;
while(getline(fin,line))
	{
	if(line[0] == '#')
		{
		continue;
		}
	if (firstLine==0)
		{
		istringstream ss1(line);
		ss1 >> N >> Na >> Nb >> Nc >> dE >> Tb >> theta;
		firstLine++;
		}
	else if (firstLine==1)
		{
		istringstream ss2(line);
		ss2 >> aq.inputChoice >> aq.fileNo >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
		ss2 >> alpha >> open >> amp;
		firstLine++;
		}
	}
fin.close();
inP = aq.inputChoice; //just because I write this a lot

//derived quantities
NT = Na + Nb + Nc;
epsilon = dE;
R = 2.0/3.0/epsilon;
alpha *= R;
vec negVec(2*N*Nb+1);
double negVal;
if (inP.compare("p") == 0)
	{
	L = 3.2*R;
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
		cout << "Tb>R so running negEig, need to have run pi with inP='b'" << endl;
		system("./negEig");
		cout << "negEig run" << endl;
		negVec = loadVector("./data/eigVec.dat",Nb);
		ifstream eigFile;
		eigFile.open("./data/eigVal.dat", ios::in);
		string lastLine = getLastLine(eigFile);
		istringstream ss(lastLine);
		double temp;
		double eigError; //should be <<1
		ss >> temp >> temp >> temp >> eigError >> negVal;
		if (eigError>1.0)
			{
			cout << "error in negEig = " << eigError << endl;
			cout << "consider restarting with different values of P and c" << endl;
			}
		}
	}
else if (inP.compare("b") == 0)
	{
	Tb = 1.5*R;
	L = 3.2*R;
	}
a = L/(N-1.0);
b = Tb/(Nb-1.0);
Ta = b*Na;
Tc = b*Nc;
Gamma = exp(-theta);

//determining number of runs
closenessA = 1.0;
closenessS = 1.0e-5;
closenessSM = 1.0e-4;
closenessD = 1.0;
closenessC = 5.0e-14;

string loop_choice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
string print_choice = aq.printChoice;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//begin loop over varying parameter
//(NB. may only be one loop)
//and defining some quantities to be used later
for (unsigned int loop=0; loop<aq.totalLoops; loop++)
	{
	//giving values of varying parameters
	if (loop_choice.compare("N") == 0)
		{
		int loopParameter = aq.minValue + (int)(aq.maxValue - aq.minValue)*loop/(aq.totalLoops-1);
		changeInt (loop_choice,loopParameter);
		}
	else if (loop_choice.compare("n")!=0)
		{
		double loopParameter = aq.minValue + (aq.maxValue - aq.minValue)*loop/(aq.totalLoops-1.0);
		changeDouble (loop_choice,loopParameter);
		}
		
	printParameters();

	//defining a time and starting the clock
	clock_t time;
	time = clock();
	
	//defining some important scalar quantities
	double S1 = 2.0/3.0; //mass of kink multiplied by lambda
	double twaction = -2.0*pi*epsilon*pow(R,2)/2.0 + 2.0*pi*R*S1;
	comp action = twaction;
	cVec erg(NT);	

	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp action_last = action;
	unsigned int runs_count = 0;
	unsigned int min_runs = 2;
	vector<double> action_test(1);	action_test[0] = 1.0;
	vector<double> sol_test(1);	sol_test[0] = 1.0;
	vector<double> solM_test(1);	solM_test[0] = 1.0;
	vector<double> delta_test(1); delta_test[0] = 1.0;
	vector<double> calc_test(1); calc_test[0] = 1.0;

	//initializing phi (=p)
	vec p(2*N*Nb+1);
	p = Eigen::VectorXd::Zero(2*N*Nb+1);
	
	//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
	double b_parameter = -1.0;
	double c_parameter = -epsilon;
	gsl_poly_solve_cubic (0, b_parameter, c_parameter, &root[0], &root[1], &root[2]);
	sort(root.begin(),root.end());

	//deterimining omega matrices for fourier transforms in spatial direction
	mat h(N,N);
	h = hFn(N,a);
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
		eigenVectors = eigensolver.eigenvectors();
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
	if (inP.compare("f")!=0)
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
				if (Tb>R)
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
		string prefix = ("./data/pi");
		string suffix = (".dat");
		string loadfile = prefix+to_string(aq.fileNo)+suffix;
		p = loadVector(loadfile,Nb);
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
	printVectorB("data/piEarly00.dat",p);
		
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
        	printVectorB("data/Chi0.dat",Chi0);
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
		erg = Eigen::VectorXcd::Zero(NT);

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
				comp tempErg =  (kineticS + pot_l + pot_e)/Dt;
				kineticS += Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0; //n.b. no contribution from time derivative term at the final time boundary
				pot_l += Dt*a*Va(Cp(j));
				pot_e += Dt*a*Vb(Cp(j));
				erg(t+Na) += + (kineticS + pot_l + pot_e)/Dt - tempErg;
				
				DDS.insert(2*j,2*j) = 1.0/b; //zero time derivative
				DDS.insert(2*j,2*(j-1)) = -1.0/b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
			else if (t==0)
				{
				comp dt = -b*i;
				comp Dt = -b*i/2.0;
				comp tempErg =  (kineticS + pot_l + pot_e)/Dt + kineticT/dt;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				kineticS += Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += Dt*a*Va(Cp(j));
				pot_e += Dt*a*Vb(Cp(j));
				erg(t+Na) += (kineticS + pot_l + pot_e)/Dt + kineticT/dt - tempErg;
				
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
				comp tempErg =  (kineticS + pot_l + pot_e)/Dt + kineticT/dt;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				kineticS += Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += Dt*a*Va(Cp(j));
				pot_e += Dt*a*Vb(Cp(j));
				erg(t+Na) += (kineticS + pot_l + pot_e)/Dt + kineticT/dt - tempErg;
				
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
                comp temp1 = a*Dt*(2.0*Cp(j)/pow(a,2.0) + dV(Cp(j)));
                comp temp2 = a*Dt*(2.0/pow(a,2.0) + ddV(Cp(j)));
                    
                minusDS(2*j) += real(temp1 - temp0*Cp(j));
                minusDS(2*j+1) += imag(temp1 - temp0*Cp(j));
                DDS.insert(2*j,2*j) = real(-temp2 + temp0);
                DDS.insert(2*j,2*j+1) = imag(temp2 - temp0);
                DDS.insert(2*j+1,2*j) = imag(-temp2 + temp0);
                DDS.insert(2*j+1,2*j+1) = real(-temp2 + temp0);
                }
            }
        action = kineticT - kineticS - pot_l - pot_e;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			if ((print_choice.compare("a")==0 || print_choice.compare("e")==0) && 1==0) //have stopped this one as it's annoying
				{
				printAction(kineticT-kineticS,pot_l,pot_e);
				}
			if (print_choice.compare("v")==0 || print_choice.compare("e")==0)
				{
				string minusDSprefix = ("./data/minusDSE");
				string minusDSsuffix = (".dat");
				string minusDSfile = minusDSprefix+inP+to_string(loop)+to_string(runs_count)+minusDSsuffix;
				printVectorB(minusDSfile,minusDS);
				}
			if (print_choice.compare("p")==0 || print_choice.compare("e")==0)
				{
				string piEarlyPrefix = ("./data/pE");
				string piEarlySuffix = (".dat");
				string piEarlyFile = piEarlyPrefix+inP+to_string(loop)+to_string(runs_count)+piEarlySuffix;
				printVectorB(piEarlyFile,p);
				}
			if (print_choice.compare("m")==0 || print_choice.compare("e")==0)
				{
				string DDSprefix = ("./data/DDSE");
				string DDSsuffix = (".dat");
				string DDSfile = DDSprefix+inP+to_string(loop)+to_string(runs_count)+DDSsuffix;
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
            -Dt0*dV(ap(l)))/dtau;
    	}
    	
    //A4.5 starting the energy and that off
    #define eta(m) real(ap(m))-root[0]
    #define ieta(m) imag(ap(m))
	cVec linErgA(Na); linErgA = Eigen::VectorXcd::Zero(Na);
	cVec linNumA(Na); linNumA = Eigen::VectorXcd::Zero(Na);

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
            -dV(ap(m));
            erg (Na-t) += a*pow(ap(m-1)-ap(m),2.0)/dtau/2.0 + dtau*pow(ap(neigh(m,1,1,Na+1))-ap(m),2.0)/a/2.0 + dtau*a*V(ap(m));
            for (unsigned int y=0; y<N; y++)
            	{
            	unsigned int l1 = t + x*(Na+1);
            	unsigned int l2 = t + y*(Na+1);
		        if (absolute(theta)<2.0e-16)
					{
		        	linErgA(Na-t) += Eomega(x,y)*eta(l1)*eta(l2) - Eomega(x,y)*ap(l1)*ap(l2);
					linNumA (Na-t) += omega(x,y)*eta(l1)*eta(l2) - omega(x,y)*ap(l1)*ap(l2);
		        	}
				else
					{
					linErgA(Na-t) += 2.0*Gamma*omega(x,y)*eta(l1)*eta(l2)/pow(1.0+Gamma,2.0) + 2*Gamma*omega(x,y)*ieta(l1)*ieta(l2)/pow(1.0-Gamma,2.0);
					linNumA(Na-t) += 2.0*Gamma*Eomega(x,y)*eta(l1)*eta(l2)/pow(1.0+Gamma,2.0) + 2.0*Gamma*Eomega(x,y)*ieta(l1)*ieta(l2)/pow(1.0-Gamma,2.0);
					}
				}
        	}
    	}
    	
    //computing boundary term, at initial time
    comp bound = 0.0;
	for (unsigned int j=0; j<N; j++)
		{
		for (unsigned int k=0; k<N; k++)
			{
			
			if (absolute(theta)<2.0e-16)
				{
				bound += 0.0;
				}
			else
				{
				unsigned int l = Na + j*(Na+1);
				unsigned int m = Na + k*(Na+1);
				bound += (1.0-Gamma)*omega(j,k)*eta(l)*eta(m)/(1.0+Gamma) + (1.0+Gamma)*omega(j,k)*ieta(l)*ieta(m)/(1.0-Gamma);
				}
			}
		}

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
            -Dt0*dV(ccp(l)))/dtau;
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
		    	erg (Na+Nb-2+t) += a*pow(ccp(l)-ccp(l-1),2.0)/dtau/2.0 + dtau*pow(ccp(neigh(l-1,1,1,Nc+1))-ccp(l-1),2.0)/a/2.0 + dtau*a*V(ccp(l-1));
		    	}
			}
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
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%16s%16s\n","runs","time","N","Na","Nb","Nc","L","Tb","R","re(action)","im(action)");
	printf("%8i%8g%8i%8i%8i%8i%8g%8g%8g%16g%16g\n",runs_count,realtime,N,Na,Nb,Nc,L,Tb,R,real(action),imag(action));
	printf("\n");
	 printf("%60s\n","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

	//printing action value
	FILE * actionfile;
	actionfile = fopen("./data/action.dat","a");
	fprintf(actionfile,"%8i%8g%8i%8i%8i%8i%8g%8g%8g%16g%16g\n",runs_count,realtime,N,Na,Nb,Nc,L,Tb,R,real(action),imag(action));
	fclose(actionfile);

	//printing output phi
	string pifile = "./data/pi"+inP+to_string(loop)+".dat";
	printVector(pifile,tp);
	gp(pifile,"repi.gp");
	
	//printing output minusDS				
	string minusDSfile = "./data/minusDS"+inP+to_string(loop)+".dat";
	printVectorB(minusDSfile,minusDS);
				
	//printing output DDS
	string DDSfile = "./data/DDS"+inP+to_string(loop)+".dat";
	printSpmat(DDSfile,DDS);
	
	//printing linErgVec
	string linErgFile = "./data/linErg"+inP+to_string(loop)+".dat";
	simplePrintCVector(linErgFile,linErgA);
	gpSimple(linErgFile);
	
	//printing erg
	string ergFile = "./data/erg"+inP+to_string(loop)+".dat";
	simplePrintCVector(ergFile,erg);
	gpSimple(ergFile);

} //closing parameter loop

return 0;
}
