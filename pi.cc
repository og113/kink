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
#include <gsl/gsl_poly.h>
#include "pf.h"

using namespace std;

int main()
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//getting variables and user inputs from inputs

aqStruct aq; //struct to hold user responses

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
		ss1 >> N >> Na >> Nb >> Nc >> dE >> theta;
		firstLine = 1;
		}
	else
		{
		istringstream ss2(line);
		ss2 >> aq.inputChoice >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
		ss2 >> alpha >> open;
		}
	}
fin.close();
inP = aq.inputChoice; //just because I write this a lot

//derived quantities
NT = Na + Nb + Nc;
epsilon = dE;
R = 2.0/3.0/epsilon;
alpha *= R;
if (inP.compare("p") == 0)
	{
	Tb = 1.0*R;
	L = 3.0*R;
	double Ltemp = 1.5*(1.5*Tb*tan(angle));
	if (Ltemp<L) //making sure to use the smaller of the two possible Ls
		{
		L=Ltemp;
		}	
	}
else if (inP.compare("b") == 0)
	{
	Tb = 1.5*R;
	L = 3.0*R;
	}
angle = asin(Tb/R);
a = L/(N-1.0);
b = Tb/(Nb-1.0);
Ta = b*(Na-1.0);
Tc = b*(Nc-1.0);

//determining number of runs
closenessA = 1.0;
closenessS = 1.0e-5;
closenessSM = 1.0e-4;
closenessD = 1.0;
closenessC = 5.0e-14;

string loop_choice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
string print_choice = aq.printChoice;

printParameters();

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

	//defining a time and starting the clock
	clock_t time;
	time = clock();
	
	//defining some important scalar quantities
	double S1 = 2.0/3.0;
	double twaction = -2.0*pi*epsilon*pow(R,2)/2.0 + 2.0*pi*R*S1;
	complex<double> action = twaction;

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
	vector<double> root(3);
	gsl_poly_solve_cubic (0, b_parameter, c_parameter, &root[0], &root[1], &root[2]);
	sort(root.begin(),root.end());

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//assigning input phi
    for (unsigned int j=0; j<N*Nb; j++)
    	{
		if (R<alpha)
			{
		    cout << "R is too small. Not possible to give thinwall input. It should be more that " << alpha;
		    }
		else
			{
			p(2*j+1) = 0.0; //imaginary parts set to zero
			if (inP.compare("b")==0)
				{
				double rho = real(sqrt(-pow(coordB(j,0),2.0) + pow(coordB(j,1),2.0))); //should be real even without real()
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
				}
			else if (inP.compare("p")==0)
				{
				double rho1 = real(sqrt(-pow(coordB(j,0),2.0) + pow(coordB(j,1)+R*cos(angle),2.0)));
				double rho2 = real(sqrt(-pow(coordB(j,0),2.0) + pow((coordB(j,1)-R*cos(angle)),2.0)));
				if ((rho1-R)<-alpha && (rho2-R)<-alpha)
					{
				    p(2*j) = root[2];
				    }
				else if ((rho1-R)>alpha || (rho2-R)>alpha)
					{
				    p(2*j) = root[0];
				    }
				else if (real(coordB(j,1))>0) //note that the coord should be real
					{
				    p(2*j) = (root[2]+root[0])/2.0 + (root[0]-root[2])*tanh((rho1-R)/2.0)/2.0;
				    }
				else if (real(coordB(j,1))<0)
					{
				    p(2*j) = (root[2]+root[0])/2.0 + (root[0]-root[2])*tanh((rho2-R)/2.0)/2.0;
				    }
				else
					{
				    p(2*j) = root[2]; //i.e. if coordB(j,1) == 0
				    }
				}
			}
		}
		
	p(2*N*Nb) = 0.5; //initializing Lagrange parameter for removing dp/dx zero mode
	
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
		vec minusDS(2*N*Nb+1);
		minusDS = Eigen::VectorXd::Zero(2*N*Nb+1); //initializing to zero
		spMat DDS(2*N*Nb+1,2*N*Nb+1);
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
		comp kinetic = 0.0;
		comp pot_l = 0.0;
		comp pot_e = 0.0;

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
				kinetic += -Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0; //n.b. no contribution from time derivative term at the final time boundary
				pot_l += Dt*a*Va(Cp(j));
				pot_e += Dt*a*Vb(Cp(j));
				
				DDS.insert(2*j,2*j) = 1.0/b; //zero time derivative
				DDS.insert(2*j,2*(j-1)) = -1.0/b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
			else if (t==0)
				{
				comp dt = -b*i;
				comp Dt = -b*i/2.0;
				kinetic += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0 - Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += Dt*a*Va(Cp(j));
				pot_e += Dt*a*Vb(Cp(j));
				
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
				kinetic += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0 - Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0;
				pot_l += Dt*a*Va(Cp(j));
				pot_e += Dt*a*Vb(Cp(j));
				
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
        action = kinetic - pot_l - pot_e;   
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			if ((print_choice.compare("a")==0 || print_choice.compare("e")==0) && 1==0) //have stopped this one as it's annoying
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
            -dV(ap(j*Na)))/dtau;
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
            -dV(ap(l));    
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
    cVec accC(N*Nc);
    accC = Eigen::VectorXcd::Zero(N*Nc);
    for (unsigned int j=0; j<N; j++)
    	{
        accC(j*Nc) = ((Dt0/pow(a,2))*(cp(neigh(j*Nc,1,1,Nc))+cp(neigh(j*Nc,1,-1,Nc))-2.0*cp(j*Nc)) \
            -dV(cp(j*Nc)))/dtau;
    	}

    //C7. run loop
    for (unsigned int j=1; j<Nc; j++)
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
		    -dV(cp(l));    
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
    vec tp(2*N*NT);
    tp = vecReal(tCp,NT*N);
    tp.conservativeResize(2*N*NT+1);
    tp(2*N*NT) = p(2*N*Nb);
    
    //stopping clock
	time = clock() - time;
	double realtime = time/1000000.0;
	
	//printing to terminal
	if (loop==0)
		{
		printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%16s%16s\n","runs","time","N","Na","Nb","Nc","L","Tb","R","re(action)","im(action)");
		}
	printf("%8i%8g%8i%8i%8i%8i%8g%8g%8g%16g%16g\n",runs_count,realtime,N,Na,Nb,Nc,L,Tb,R,real(action),imag(action));

	//printing action value
	FILE * actionfile;
	actionfile = fopen("./data/action.dat","a");
	fprintf(actionfile,"%8i%8g%8i%8i%8i%8i%8g%8g%8g%16g%16g\n",runs_count,realtime,N,Na,Nb,Nc,L,Tb,R,real(action),imag(action));
	fclose(actionfile);

	//printing output phi
	string oprefix = ("./data/pi");
	string osuffix = (".dat");
	string outfile = oprefix+to_string(loop)+osuffix;
	printVector(outfile,tp);

} //closing parameter loop

return 0;
}
