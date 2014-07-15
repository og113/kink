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
printParameters();

aqStruct aq{ //struct to hold user responses (answers to questions)
aq.fileNo;
aq.maxTheta;
aq.totalLoops = 1;
aq.loopChoice = 'n';
aq.minValue;
aq.maxValue;
aq.printChoice = "n";
aq.printRun = -1;
};

askQuestions(aq);

loopChoice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
printChoice = aq.printChoice;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//begin loop over varying parameter
//(NB. may only be one loop)
//and defining some quantities to be used later
for (int loop=0; loop<aq.totalLoops; loop++)
	{
	//giving values of varying parameters
	if (loopChoice.compare(0,1,"N") == 0)
		{
		int loopParameter = aq.minValue + (int)(aq.maxValue - aq.minValue)*loop/(aq.totalLoops-1);
		changeInt (loopChoice,loopParameter);
		}
	else if (loopChoice.compare("n")!=0)
		{
		double loopParameter = min_value + (max_value - min_value)*loop/(total_loops-1.0);
		changeDouble (loopChoice,loopParameter);
		}

	//defining a time
	clock_t time;
	clock_t wait;
	
	//defining some important scalar quantities
	complex<double> action = 2.0;
	double S_1 = 2.0*pow(mass,3)/3.0/lambda;
	double twaction = -2.0*pi*epsilon*pow(R,2)/2.0 + 2.0*pi*R*S_1;
	int alpha = 5; //gives span over which tanh is used

	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp actionLast = 1.0;
	int runs_count = 0;
	double action_test = 1.0;
	double vector_test = 1.0;
	unsigned int min_runs = 3;

	//initializing phi (=p)
	vec p(2*N*Nb+1);
	
	//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
	double b_parameter = -pow(v,2.0);
	double c_parameter = epsilon/v/lambda;
	vector<double> root(3);
	gsl_poly_solve_cubic (0, b_parameter, c_parameter, &root[0], &root[1], &root[2]);
	sort(root.begin(),root.end());

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//assigning input phi
    for (unsigned int j=0, j<N*Nb; j++)
    	{
		complex<double> rho1Sqrd = -pow(coordB(j,0),2) + pow(coordB(j,1)+R*cos(angle)),2);
		complex<double> rho2Sqrd = -pow(coordB(j,0),2) + pow((coordB(j,1)-R*cos(angle)),2); 
		double rho1 = sqrt(real(rho1Sqrd)); //should be real even without real()
		double rho2 = sqrt(real(rho2Sqrd));
		if (R<alpha/mass)
			{
		    cout << "X = R*mass is too small. not possible to give thinwall input. it should be less that " << alpha;
		    }
		else
			{
		    if (rho1<(R-alpha/mass) && rho2<(R-alpha/mass))
		    	{
		        p(2*j+1) = root(1);
		        }
		    else if (rho1>(R+alpha/mass) || rho2>(R+alpha/mass))
		    	{
		        p(2*j+1) = root(3);
		        }
		    else if (real(coordB(j,1))>0) //note that the coord should be real
		    	{
		        p(2*j+1) = v*tanh(mass*(rho1-R)/2);
		        }
		    else if (real(eCoord(j,1))<0)
		    	{
		        p(2*j+1) = v*tanh(mass*(rho2-R)/2);
		        }
		    else
		    	{
		        p(2*j+1) = root(1); //i.e. if eCoord(j,1) == 0
		        }
			}
		}
		
	p(N*Nb+1) = v; //initializing Lagrange parameter for removing dp/dx zero mode
	
	//fixing input periodic instanton to have zero time derivative at time boundaries
    double open = 0.5;//value of 0 assigns all weight to boundary, value of 1 to neighbour of boundary
    for (unsigned int j=0;j<N;j++)
    	{
        p(2*j*Nb+1) = open*p(2*j*Nb+1) + (1-open)*p(2*(j*Nb+1)+1);%intiial time real
        p(2*(j*Nb+1)+1) = open*p(2*j*Nb+1) + (1-open)*p(2*(j*Nb+1)+1);
        p(2*j*Nb+2) = open*p(2*j*Nb+2) + (1-open)*p(2*(j*Nb+1)+2);%initial time imag
        p(2*(j*Nb+1)+2) = open*p(2*j*Nb+2) + (1-open)*p(2*(j*Nb+1)+2);
        p(2*((j+1)*Nb-1)) = open*p(2*((j+1)*Nb)) + (1-open)*p(2*((j+1)*Nb-1));%final time real
        p(2*((j+1)*Nb)) = open*p(2*((j+1)*Nb)) + (1-open)*p(2*((j+1)*Nb-1));
        p(2*((j+1)*Nb-2)+2) = open*p(2*((j+1)*Nb-1)+2) + (1-open)*p(2*((j+1)*Nb-2)+2);%final time imag
        p(2*((j+1)*Nb-1)+2) = open*p(2*((j+1)*Nb-1)+2) + (1-open)*p(2*((j+1)*Nb-2)+2);
		}
		
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
		runs_count ++;
		action_last = action;

		// allocating memory for DS, DDS
		vec minusDS(2*N*Nb+1);
		minusDS = Eigen::VectorXd::Zero(2*N*Nb+1);
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
		
		//defining the zero mode at the final time boundary
		cVec Chi0(N);
		for (unsigned int j=0; j<(N-1); j++)
			{
            Chi0(j) = p(2*((j+2)*Nb-1))+1i*p(2*((j+2)*Nb-1))-p(2*((j+1)*Nb-1))-1i*p(2*((j+1)*Nb-1)); //note only use real derivative - this is a fudge due to initial input
             }
        Chi0(N-1) = p(2*(Nb-1))+1i*p(2*(Nb-1))-p(2*(N*Nb-1))-1i*p(2*(N*Nb-1)); //written directly to avoid using neigh
        norm = Chi0.dot(Chi0);
        Chi0 = Chi0/norm;
		
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
			unsigned int x = intCoord(j,1,Nb);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//boundaries			
			if (t==(Nb-1))
				{
				comp Dt = -b*i/2.0;
				kinetic += -Dt*pow(Cp(neigh(j,k,1,Nb))-Cp(j),2)/a/2.0; //n.b. no contribution from time derivative term at the final time boundary
				pot_l += -Dt*a*lambda*pow(pow(Cp(j),2)-pow(v,2),2)/8.0;
				pot_e += -Dt*a*epsilon*(Cp(j)-v)/v/2.0;
				
				DDS.insert(2*j,2*j) = 1.0/b; //zero time derivative
				DDS.insert(2*j,2*(j-1)) = -1.0/b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				
				DDS.insert(2*j+1,2*N*Nb+1) = real(siteMeasure*Chi0(x+1)); //zero mode lagrange constraint
                DDSm(2*j+2,2*N*Nb+1) = imag(siteMeasure*Chi0(x+1)); //the constraint is real but its derivative wrt phi may be complex
				}
			else if (t==0)
				{
				comp dt = -b*i;
				comp Dt = -b*i/2.0;
				kinetic += a*pow(Cp(j+1)-Cp(j),2)/dt/2.0\
				-Dt*pow(Cp(neigh(j,k,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += -Dt*a*lambda*pow(pow(Cp(j),2)-pow(v,2),2)/8.0;
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
				kinetic += a*pow(Cp(j+1)-Cp(j),2)/dt/2.0\
				-Dt*pow(Cp(neigh(j,k,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += -Dt*a*lambda*pow(pow(Cp(j),2)-pow(v,2),2)/8.0;
				pot_e += -Dt*a*epsilon*(Cp(j)-v)/v/2.0;
				
                for (unsigned int k=0; k<2*2; k++)
                	{
                    sign = (-1)^k;
                    direc = floor(k/2);
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
                        neighb = neigh(j,direc,sign,Nb);
                        minusDS(2*j) += - real(Dtj*Cp(neighb)/a);
                        minusDS(2*j+1) += - imag(Dtj*Cp(neighb)/a);
                        DDS.insert(2*j,2*neighb) = real(Dt/a);
                        DDS.insert(2*j,2*neighb+1) = -imag(Dt/a);
                        DDS.insert(2*j+1,2*neighb) = imag(Dt/a);
                        DDS.insert(2*j+1,2*neighb+1) = real(Dt/a);
                        }
                    }
                comp temp0 = 2*a*dt;
                comp temp1 = a*Dt*(2*Cp(j)/pow(a,2) + (lambda/2.0)*Cp(j)*(pow(Cp(j),2)-pow(v,2)) + epsilon/2.0/v);
                comp temp2 = a*Dt*(2/pow(a,2) + (lambda/2.0)*(3.0*pow(Cp(j),2) - pow(v,2)));
                    
                minusDS(2*j+1) += real(temp1 - temp0*Cp(j+1));
                minusDS(2*j+2) += imag(temp1 - temp0*Cp(j+1));
                DDS.insert(2*j,2*j) = real(-temp2 + temp0);
                DDS.insert(2*j,2*j+1) = imag(temp2 - temp0);
                DDS.insert(2*j+1,2*j) = imag(-temp2 + temp0);
                DDS.insert(2*j+1,2*j+1) = real(-temp2 + temp0);
                }
            }

        for (unsigned int j=0; j<N; j++) //lagrange multiplier terms
        	{
            minusDS(2*N*Nb) = minusDS(2*N*Nb) - real(a*b*Chi0(j)*Cp((j+1)*Nb-1));
            DDS.insert(2*N*Nb,2*((j+1)*Nb-1)) = real(a*b*Chi0(j+1));
            DDS.insert(2*N*Nb,2*((j+1)*Nb-1)+1) = -imag(a*b*Chi0(j+1));
            }
        action = kinetic + potL + potE;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//printing early if desired
		if (runs_count == aq.printRun)
			{
			if (printChoice.compare("a"))
				{
				printAction(kinetic,potL,potE)
				}
			else if (printChoice.compare("v"))
				{
				printVector("data/minusDS.dat",minusDS);
				}
			else if (printChoice.compare("p"))
				{
				printVector("data/piEarly.dat",p);
				}
			else if (printChoice.compare("m"))
				{
				printSpmat("data/DDS.dat",DDS);
				}
			}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//solving for delta using the Newton-Raphson method
		vec delta(2*N*Nb+1);
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
		
		char stop_wait;
		char print_wait;
		clock_t stop_time = clock();
		bool bool_wait = false; //set to false if you want program to stop if the looping is slow and ask the user whether or not they want to print
		convergence(runs_count,action_test,vector_test,time,&wait,&print_wait,&print_choice, &print_run, kinetic, pot_l, pot_e, bool_wait);
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
    dtau = -b;
    Dt0 = dtau; //b/2*(-1+1i*up); - this is surely wrong!!
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
            l = j+k*Na;
            velA(l) = velA(l-1) + dtau*accA(l-1);
            ap(l) = ap(l-1) + dtau*velA(l);
        	}
        for (unsigned int k=0; k<N; k++)
        	{
            l = j+k*Na;
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
		    l = j+k*Nc;
		    vel(l) = vel(l-1) + dtau*acc(l-1);
		    ap(l) = ap(l-1) + dtau*vel(l);
			}
		for (unsigned int k=0; k<N; k++)
			{
		    l = j+k*Nc;
		    acc(l) = (1.0/pow(a,2))*(cp(neigh(l,1,1,Nc))+cp(neigh(l,1,-1,Nc))-2.0*cp(l)) \
		    -(lambda/2.0)*ap(l)*(pow(cp(l),2)-pow(v,2)) - epsilon/2.0/v;    
			}
		}
    

    //12. combine phi with ap and cp and save combination to file
    cVec tCp(NT*N);
    for (unsigned int j=0; j<NT*N; j++)
    	{
        t = intCoord(j,0,NT);
        x = intCoord(j,1,NT);
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
    vec t(N*NT)
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
	printVec(outfile,tp);

} //closing parameter loop

return 0;
}
