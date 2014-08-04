//program to generate the periodic instantons
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
#include <gsl/gsl_poly.h>
#include "pf.h"

using namespace std;

int main()
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//user interface
if ( L > Ltemp ) //making sure to use the smaller of the two possible Ls
	{
	L = Ltemp;
	a = L/(N-1.0);
	}
printParameters();

aqStruct aq; //struct to hold user responses (answers to questions)
aq.fileNo = 0;
aq.maxTheta = 0;
aq.totalLoops = 1;
aq.loopChoice = "n";
aq.minValue = 0;
aq.maxValue = 0;
aq.printChoice = "n";
aq.printRun = -1;

askQuestions(aq);

string loop_choice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
string print_choice = aq.printChoice;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//begin loop over varying parameter
//(NB. may only be one loop)
//and defining some quantities to be used later
for (unsigned int loop=0; loop<aq.totalLoops; loop++)
	{
	//giving values of varying parameters
	if (loop_choice.compare(0,1,"N") == 0)
		{
		int loopParameter = aq.minValue + (int)(aq.maxValue - aq.minValue)*loop/(aq.totalLoops-1);
		changeInt (loop_choice,loopParameter);
		}
	else if (loop_choice.compare("n")!=0)
		{
		double loopParameter = aq.minValue + (aq.maxValue - aq.minValue)*loop/(aq.totalLoops-1.0);
		changeDouble (loop_choice,loopParameter);
		}

	//defining a time
	clock_t time;
	clock_t wait;
	
	//starting the clock
	time = clock();
	wait = clock();
	
	//defining some important scalar quantities
	complex<double> action = 2.0;
	//double S_1 = 2.0*pow(mass,3)/3.0/lambda;
	//double twaction = -2.0*pi*epsilon*pow(R,2)/2.0 + 2.0*pi*R*S_1;
	double alpha = 8.0; //gives span over which tanh is used

	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp action_last = 1.0;
	unsigned int runs_count = 0;
	double action_test = 1.0;
	double vector_test = 1.0;
	unsigned int min_runs = 3;

	//initializing phi (=p)
	vec p(2*N*Nb+1);
	p = Eigen::VectorXd::Zero(2*N*Nb+1);
	
	//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
	double b_parameter = -pow(v,2.0);
	double c_parameter = epsilon/v/lambda;
	vector<double> root(3);
	gsl_poly_solve_cubic (0, b_parameter, c_parameter, &root[0], &root[1], &root[2]);
	sort(root.begin(),root.end());

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//assigning input phi
    for (unsigned int j=0; j<N*Nb; j++)
    	{
		complex<double> rho1Sqrd = -pow(coordB(j,0),2.0) + pow(coordB(j,1)+R*cos(angle),2.0);
		complex<double> rho2Sqrd = -pow(coordB(j,0),2.0) + pow((coordB(j,1)-R*cos(angle)),2.0); 
		double rho1 = sqrt(real(rho1Sqrd)); //should be real even without real()
		double rho2 = sqrt(real(rho2Sqrd));
		if (R<alpha/mass)
			{
		    cout << "X = R*mass is too small. not possible to give thinwall input. it should be less that " << alpha;
		    }
		else
			{
			p(2*j+1) = 0.0; //imaginary parts set to zero
		    if (rho1<(R-alpha/mass) && rho2<(R-alpha/mass))
		    	{
		        p(2*j) = root[0];
		        }
		    else if (rho1>(R+alpha/mass) || rho2>(R+alpha/mass))
		    	{
		        p(2*j) = root[2];
		        }
		    else if (real(coordB(j,1))>0) //note that the coord should be real
		    	{
		        p(2*j) = (root[2]+root[0])/2.0 + (root[2]-root[0])*tanh(mass*(rho1-R)/2.0)/2.0;
		        }
		    else if (real(coordB(j,1))<0)
		    	{
		        p(2*j) = (root[2]+root[0])/2.0 + (root[2]-root[0])*tanh(mass*(rho2-R)/2.0)/2.0;
		        }
		    else
		    	{
		        p(2*j) = root[0]; //i.e. if coordB(j,1) == 0
		        }
			}
		}
		
	p(N*Nb) = v; //initializing Lagrange parameter for removing dp/dx zero mode
	
	//fixing input periodic instanton to have zero time derivative at time boundaries
    double open = 1.0;//value of 0 assigns all weight to boundary, value of 1 to neighbour of boundary
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
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//beginning newton-raphson loop
	bool Xwait = true;
	while ((action_test > closenessA || runs_count<min_runs || vector_test>closenessV) && Xwait )
		{
		//quantities used to stop newton-raphson loop
		action_test = abs(action - action_last)/abs(action_last);
		if (action_test < closenessA)
			{
			break;
			}
		runs_count ++;
		action_last = action;

		// allocating memory for DS, DDS
		vec minusDS(2*N*Nb+1);
		minusDS = Eigen::VectorXd::Zero(2*N*Nb+1); //initializing to zero
		spMat DDS(2*N*Nb+1,2*N*Nb+1);
		Eigen::VectorXi DDS_to_reserve(2*N*Nb+1);//number of non-zero elements per column
		DDS_to_reserve(0) = 2;
		DDS_to_reserve(1) = 2;
		DDS_to_reserve(2*N*Nb-2) = 2;
		DDS_to_reserve(2*N*Nb-1) = 2;
		DDS_to_reserve(2*N*Nb) = 2*N;
		for (lint j=1;j<(N*Nb-1);j++)
			{
			DDS_to_reserve(2*j) = 2*(2*2+1);
			DDS_to_reserve(2*j+1) = 2*(2*2+1);
			}
		DDS.reserve(DDS_to_reserve);
		
		//defining the zero mode at the final time boundary and the time step before
		vec Chi0(Nb*N);
		Chi0 = Eigen::VectorXd::Zero(N*Nb);
		for (unsigned int j=0; j<N; j++)
			{
			unsigned int pos = (j+1)*Nb-1;
            Chi0(pos) = p(2*neigh(pos,1,1,Nb))-p(2*neigh(pos,1,-1,Nb)); //final time slice
            //Chi0(pos-1) = p(2*neigh(pos-1,1,1,Nb))-p(2*neigh(pos-1,1,-1,Nb)); //penultimate time slice
            }
        //double norm; //normalizing
        //norm = Chi0.dot(Chi0);
        //norm = pow(norm,0.5);
        //Chi0 /= norm;
        //Chi0 *= v;
		
		//initializing to zero
		comp kinetic = 0.0;
		comp pot_l = 0.0;
		comp pot_e = 0.0;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//assigning values to minusDS and DDS and evaluating action
		for (unsigned long int j = 0; j < N*Nb; j++)
			{		
			unsigned int t = intCoord(j,0,Nb); //coordinates
			//unsigned int x = intCoord(j,1,Nb); //currently unused
			
			if (Chi0(j)>1.0e-16) //zero mode lagrange constraint
				{
				DDS.insert(2*j,2*N*Nb) = a*Chi0(j); 
				DDS.insert(2*N*Nb,2*j) = a*Chi0(j);
		    	minusDS(2*j) += -a*Chi0(j)*p(2*N*Nb);
		    	minusDS(2*N*Nb) +=- a*Chi0(j)*p(2*j-1);
		    	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//boundaries			
			if (t==(Nb-1))
				{
				comp Dt = -b*i/2.0;
				kinetic += -Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0; //n.b. no contribution from time derivative term at the final time boundary
				pot_l += -Dt*a*lambda*pow(pow(Cp(j),2.0)-pow(v,2.0),2.0)/8.0;
				pot_e += -Dt*a*epsilon*(Cp(j)-v)/v/2.0;
				
				DDS.insert(2*j,2*j) = 1.0/b; //zero time derivative
				DDS.insert(2*j,2*(j-1)) = -1.0/b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
			else if (t==0)
				{
				comp dt = -b*i;
				comp Dt = -b*i/2.0;
				kinetic += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0\
				-Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += -Dt*a*lambda*pow(pow(Cp(j),2.0)-pow(v,2.0),2.0)/8.0;
				pot_e += -Dt*a*epsilon*(Cp(j)-v)/v/2.0;
				
				DDS.insert(2*j,2*j) = -1.0/b; //zero time derivative
				DDS.insert(2*j,2*(j+1)) = 1.0/b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//bulk
			else
				{
				comp dt = -b*i;
				comp Dt = -b*i;
				kinetic += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0\
				-Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += -Dt*a*lambda*pow(pow(Cp(j),2.0)-pow(v,2.0),2.0)/8.0;
				pot_e += -Dt*a*epsilon*(Cp(j)-v)/v/2.0;
				
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
                        int neighb = neigh(j,direc,sign,Nb);
                        
                        minusDS(2*j) += - real(Dt*Cp(neighb)/a);
                        minusDS(2*j+1) += - imag(Dt*Cp(neighb)/a);
                        DDS.insert(2*j,2*neighb) = real(Dt/a);
                        DDS.insert(2*j,2*neighb+1) = -imag(Dt/a);
                        DDS.insert(2*j+1,2*neighb) = imag(Dt/a);
                        DDS.insert(2*j+1,2*neighb+1) = real(Dt/a);
                        }
                    }
                comp temp0 = 2.0*a/dt;
                comp temp1 = a*Dt*(2.0*Cp(j)/pow(a,2.0) + (lambda/2.0)*Cp(j)*(pow(Cp(j),2.0)-pow(v,2.0)) + epsilon/2.0/v);
                comp temp2 = a*Dt*(2.0/pow(a,2.0) + (lambda/2.0)*(3.0*pow(Cp(j),2.0) - pow(v,2.0)));
                    
                minusDS(2*j) += real(temp1 - temp0*Cp(j));
                minusDS(2*j+1) += imag(temp1 - temp0*Cp(j));
                DDS.insert(2*j,2*j) = real(-temp2 + temp0);
                DDS.insert(2*j,2*j+1) = imag(temp2 - temp0);
                DDS.insert(2*j+1,2*j) = imag(-temp2 + temp0);
                DDS.insert(2*j+1,2*j+1) = real(-temp2 + temp0);
                }
            }
        action = kinetic + pot_l + pot_l;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			if (print_choice.compare("a")==0 || print_choice.compare("e")==0)
				{
				printAction(kinetic,pot_l,pot_e);
				}
			if (print_choice.compare("v")==0 || print_choice.compare("e")==0)
				{
				string minusDSprefix = ("./data/minusDS");
				string minusDSsuffix = (".dat");
				string minusDSfile = minusDSprefix+to_string(loop)+to_string(runs_count)+minusDSsuffix;
				printVectorB(minusDSfile,minusDS);
				}
			if (print_choice.compare("p")==0 || print_choice.compare("e")==0)
				{
				string piEarlyPrefix = ("./data/piEarly");
				string piEarlySuffix = (".dat");
				string piEarlyFile = piEarlyPrefix+to_string(loop)+to_string(runs_count)+piEarlySuffix;
				printVectorB(piEarlyFile,p);
				}
			if (print_choice.compare("m")==0 || print_choice.compare("e")==0)
				{
				string DDSprefix = ("./data/DDS");
				string DDSsuffix = (".dat");
				string DDSfile = DDSprefix+to_string(loop)+to_string(runs_count)+DDSsuffix;
				printSpmat(DDSfile,DDS);
				}
			}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//solving for delta using the Newton-Raphson method
		double small = minusDS.dot(minusDS);
		small = pow(small,0.5);
		vector_test = small/(2*N*Nb+1); //gives average size of elements of minusDS
		if (vector_test < closenessV)
			{
			break; //solution has been found
			}
		
		vec delta(2*N*Nb+1);
		delta = Eigen::VectorXd::Zero(2*N*Nb+1);
		DDS.makeCompressed();
		Eigen::SparseLU<spMat> solver;
		
		solver.analyzePattern(DDS);
		if(solver.info()!=Eigen::Success)
			{
			cout << "DDS pattern analysis failed" << endl;
			return 0;
			}
			
		solver.factorize(DDS);
		if(solver.info()!=Eigen::Success) 
			{
			cout << "LU factorization failed" << endl;
			return 0;
			}
		delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side
		if(solver.info()!=Eigen::Success)
			{
			cout << "solving failed" << endl;
			return 0;
			}

		//assigning values to phi
		p += delta;
		
		//passing changes on to complex vector
		Cp = vecComplex(p,N*Nb);
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//convergence issues
		vector_test = minusDS.maxCoeff();
		vector_test = absolute(vector_test);
		
		char print_wait;
		bool bool_wait = false; //set to false if you want program to stop if the looping is slow and ask the user whether or not they want to print
		convergence(runs_count,action_test,vector_test,time,&wait,&print_wait,&print_choice, &aq.printRun, kinetic, pot_l, pot_e, bool_wait);
		} //closing "runs" while loop
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//propagating solution along minkowskian time paths
	
	//propagating solution back in minkowskian time
    //A1. initialize mp==mphi using last point of ephi and zeros- use complex phi
    cVec ap(N*Na); //phi on section "a"
    ap = Eigen::VectorXcd::Zero(N*Na);
    for (unsigned int j=0; j<N; j++)
    	{
        ap(j*Na) = Cp(j*Nb);
    	}

    //A2. initialize vel - defined at half steps, first step being at t=-1/2,
    //vel(t+1/2) := (p(t+1)-p(t))/dt
    cVec velA (N*Na);
    velA = Eigen::VectorXcd::Zero(N*Na);
    double dtau = -b;
    double Dt0 = dtau; //b/2*(-1+1i*up); - this is surely wrong!!
    for (unsigned int j=0; j<N; j++)
    	{
        velA(j*Na) = 0; //due to boundary condition
    	}

    
    //A3. initialize acc using phi and expression from equation of motion and zeros-complex
    cVec accA(N*Na);
    accA = Eigen::VectorXcd::Zero(N*Na);
    for (unsigned int j=0; j<N; j++)
    	{
        accA(j*Na) = ((Dt0/pow(a,2))*(ap(neigh(j*Na,1,1,Na))+ap(neigh(j*Na,1,-1,Na))-2.0*ap(j*Na)) \
            -(lambda*Dt0/2.0)*ap(j*Na)*(pow(ap(j*Na),2)-pow(v,2)) - epsilon*Dt0/2/v)/dtau;
    	}

    //A7. run loop
    for (unsigned int j=1; j<Na; j++)
    	{
        for (unsigned int k=0; k<N; k++)
        	{
            unsigned int l = j+k*Na;
            velA(l) = velA(l-1) + dtau*accA(l-1);
            ap(l) = ap(l-1) + dtau*velA(l);
        	}
        for (unsigned int k=0; k<N; k++)
        	{
            unsigned int l = j+k*Na;
            accA(l) = (1.0/pow(a,2))*(ap(neigh(l,1,1,Na))+ap(neigh(l,1,-1,Na))-2.0*ap(l)) \
            -(lambda/2.0)*ap(l)*(pow(ap(l),2)-pow(v,2)) - epsilon/2.0/v;    
        	}
    	}

    //now propagating forwards along c
    //C2. initialize mp==mphi using last point of ephi and zeros- use complex phi
    cVec cp(N*Nc); //phi on section "c"
    cp = Eigen::VectorXcd::Zero(N*Nc);
    for (unsigned int j=0; j<N; j++)
    	{
        cp(j*Nc) = Cp(j*Nb+Nb-1);
    	}

    //C3. initialize vel - defined at half steps, first step being at t=-1/2,
    //vel(t+1/2) := (p(t+1)-p(t))/dt
    cVec velC (N*Nc);
    velC = Eigen::VectorXcd::Zero(N*Nc);
    dtau = b;
    Dt0 = dtau; //b/2*(-1+1i*up); - this is surely wrong!!
    for (unsigned int j=0; j<N; j++)
    	{
        velC(j*Nc) = 0; //due to boundary condition
    	}

    //C4. initialize acc using phi and expression from equation of motion and zeros-complex
    cVec accC(N*Na);
    accC = Eigen::VectorXcd::Zero(N*Na);
    for (unsigned int j=0; j<N; j++)
    	{
        accC(j*Nc) = ((Dt0/pow(a,2))*(cp(neigh(j*Nc,1,1,Nc))+cp(neigh(j*Nc,1,-1,Nc))-2.0*cp(j*Nc)) \
            -(lambda*Dt0/2.0)*cp(j*Na)*(pow(cp(j*Na),2)-pow(v,2)) - epsilon*Dt0/2/v)/dtau;
    	}

    //C7. run loop
    for (unsigned int j=1; j<Na; j++)
		{
		for (unsigned int k=0; k<N; k++)
			{
		    unsigned int l = j+k*Nc;
		    velC(l) = velC(l-1) + dtau*accC(l-1);
		    ap(l) = ap(l-1) + dtau*velC(l);
			}
		for (unsigned int k=0; k<N; k++)
			{
		    unsigned int l = j+k*Nc;
		    accC(l) = (1.0/pow(a,2))*(cp(neigh(l,1,1,Nc))+cp(neigh(l,1,-1,Nc))-2.0*cp(l)) \
		    -(lambda/2.0)*ap(l)*(pow(cp(l),2)-pow(v,2)) - epsilon/2.0/v;    
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
            t = Na-1-t;
            tCp(j) = ap(t+x*Na);
            }
        else if (t<(Na+Nb))
        	{
            t = t - Na;
            tCp(j) = Cp(t+x*Nb);
            }
        else
        	{
            t = t - Na - Nb;
            tCp(j) = cp(t+x*Nc);
        	}
    	}
    	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	    	//misc end of program tasks - mostly printing
    
    //making real vec from complex one
    vec tp(N*NT);
    tp = vecReal(tCp,NT*N);
    tp.conservativeResize(N*NT+1);
    tp(2*N*NT) = p(2*N*Nb-1);
    
    //stopping clock
	time = clock() - time;
	double realtime = time/1000000.0;
	
	//printing to terminal
	if (loop==0)
		{
		printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%16s%16s\n","runs","time","N","Na","Nb","Nc","L","Tb","R","mass","lambda","re(action)","im(action)");
		}
	printf("%8i%8g%8i%8i%8i%8i%8g%8g%8g%8g%8g%16g%16g\n",runs_count,realtime,N,Na,Nb,Nc,L,Tb,R,mass,lambda,real(action),imag(action));

	//printing action value
	FILE * actionfile;
	actionfile = fopen("./data/action.dat","a");
	fprintf(actionfile,"%8i%8g%8i%8i%8i%8i%8g%8g%8g%8g%8g%16g%16g\n",runs_count,realtime,N,Na,Nb,Nc,L,Tb,R,mass,lambda,real(action),imag(action));
	fclose(actionfile);

	//printing output phi
	string oprefix = ("./data/pi");
	string osuffix = (".dat");
	string outfile = oprefix+to_string(loop)+osuffix;
	printVector(outfile,tp);

} //closing parameter loop

return 0;
}
