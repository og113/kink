#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "pf.h"

int main()
{
system("dir ./data/* > dataFiles");
ifstream fData;
fData.open("dataFiles", ios::in);
string line;
vector<string> filenames;
while(getline(fData,line))
	{
	if(line[0] == '#')
		{
		continue;
		}
	if (!line.empty())
		{
		string temp1, temp2;
		istringstream ss(line);
		ss >> temp1 >> temp2;
		filenames.push_back(temp1);
		filenames.push_back(temp2);
		}
	}
fData.close();

ifstream fmainin;
fmainin.open("mainInputs", ios::in);
lint minFile;
lint maxFile;
while(getline(fmainin,line))
	{
	if(line[0] == '#')
		{
		continue;
		}
	if (!line.empty())
		{
		istringstream ss(line);
		ss >> minFile >> maxFile;
		}
	}
fmainin.close();


aqStruct aq; //struct to hold user responses

//taking in information from inputs
ifstream fin;
fin.open("inputs", ios::in);
unsigned int lineInt = 0;
while(getline(fin,line))
	{
	if(line[0] == '#')
		{
		continue;
		}
	if (lineInt==0 && !line.empty())
		{
		istringstream ss1(line);
		ss1 >> N >> Na >> Nb >> Nc >> dE >> Tb >> theta;
		lineInt++;
		}
	else if (!line.empty())
		{
		istringstream ss2(line);
		ss2 >> aq.inputChoice >> aq.fileNo >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
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
	L = 3.0*R;
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
		cout << "Tb>R, cannot define angle" << endl;
		}
	}
else if (inP.compare("b") == 0)
	{
	Tb = 1.5*R;
	L = 3.0*R;
	}
a = L/(N-1.0);
b = Tb/(Nb-1.0);
Ta = b*Na;
Tc = b*Nc;

//determining number of runs
closenessA = 1.0;
closenessS = 1.0e-5;
closenessSM = 1.0e-4;
closenessD = 1.0;
closenessC = 5.0e-14;

string loop_choice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
string print_choice = aq.printChoice;

//loading initial phi
vec p(2*NT+2);
string inPrefix = ("./data/pi");
string inSuffix = (".dat");
string inFile = inPrefix + to_string(aq.fileNo) + inSuffix;
p = loadVector(inFile,NT);

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
	
	//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
	double b_parameter = -1.0;
	double c_parameter = -epsilon;
	vector<double> root(3);
	gsl_poly_solve_cubic (0, b_parameter, c_parameter, &root[0], &root[1], &root[2]);
	sort(root.begin(),root.end());
	
	//very early vector print
	printVector("data/mainEarly00.dat",p);
		
	//defining complexified vector Cp
	cVec Cp(NT*N);
	Cp = vecComplex(p,N*NT);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//beginning newton-raphson loop
	while ((sol_test.back()>closenessS || solM_test.back()>closenessSM || runs_count<min_runs))
		{
		runs_count ++;
		
		//defining the zero mode at the final time boundary and the time step before
		vec chiX(NT*N);	chiX = Eigen::VectorXd::Zero(N*NT); //to fix spatial zero mode
		vec chiT(NT*N);	chiT = Eigen::VectorXd::Zero(N*NT); //to fix real time zero mode
		for (unsigned int j=0; j<N; j++)
			{
			unsigned int posC = j*NT+(Na+Nb-1);
			unsigned int posA = j*NT;
			unsigned int posDm = j*NT+(NT-2);
            chiX(posC) = p(2*neigh(posC,1,1,NT))-p(2*neigh(posC,1,-1,NT));
            chiT(posA) = p(2*(posA+1))-p(2*posA);
            chiT(posD) = p(2*(posDm+1))-p(2*posDm);
            }
        if (runs_count==1) //printing chis
        	{
        	printVector("data/chiX.dat",chiX);
        	printVector("data/chiT.dat",chiT);
        	}
        	
        // allocating memory for DS, DDS
		vec minusDS(2*N*NT+1);
		minusDS = Eigen::VectorXd::Zero(2*N*NT+1); //initializing to zero
		spMat DDS(2*N*NT+1,2*N*NT+1);
		DDS.setZero(); //just making sure
		Eigen::VectorXi DDS_to_reserve(2*N*NT+1);//number of non-zero elements per column
		DDS_to_reserve = Eigen::VectorXi::Constant(2*N*NT+1,11);
		DDS_to_reserve(0) = 3; //these need to be changed when boundary conditions need to be more compicated
		DDS_to_reserve(1) = 3;
		DDS_to_reserve(2*N*NT-2) = 3;
		DDS_to_reserve(2*N*NT-1) = 3;
		DDS_to_reserve(2*N*NT) = N;
		DDS.reserve(DDS_to_reserve);
		
		//initializing to zero
		comp kineticS = 0.0;
		comp kineticT = 0.0;
		comp pot_l = 0.0;
		comp pot_e = 0.0;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//assigning values to minusDS and DDS and evaluating action
		for (unsigned long int j = 0; j < N*NT; j++)
			{		
			unsigned int t = intCoord(j,0,NT); //coordinates
			comp Dt = DtFn(t);
			comp dt = dtFn(t);
			
			if (absolute(chiX(j))>2.0e-16) //spatial zero mode lagrange constraint
				{
				DDS.insert(2*j,2*N*NT) = a*chiX(j); 
				DDS.insert(2*N*NT,2*j) = a*chiX(j);
		    	minusDS(2*j) += -a*chiX(j)*p(2*N*NT);
		    	minusDS(2*N*NT) += -a*chiX(j)*p(2*j);
		    	}
		    	
		    if (absolute(chiT(j))>2.0e-16)
		    	{
		    	DDS.insert(2*(j+1),2*N*NT) = a*chiT(j); 
		    	DDSm(c3) = 2*(j+1)+1; DDSn(c3) = 2*Tdim+2; DDSv(c3) = a*chiT(2*j+1); %there should be no chiT on the final time slice or this line will go wrong
                DDSm(c3) = 2*Tdim+2; DDSn(c3) = 2*(j+1)+1; DDSv(c3) = a*chiT(2*j+1);
                DDSm(c3) = 2*j+1; DDSn(c3) = 2*Tdim+2; DDSv(c3) = -a*chiT(2*j+1);
                DDSm(c3) = 2*Tdim+2; DDSn(c3) = 2*j+1; DDSv(c3) = -a*chiT(2*j+1);
                minusDS(2*(j+1)+1) = minusDS(2*(j+1)+1) - a*chiT(2*j+1)*p(2*Tdim+2);
                minusDS(2*Tdim+2) = minusDS(2*Tdim+2) - a*chiT(2*j+1)*p(2*j+1);
                minusDS(2*j+1) = minusDS(2*j+1) + a*chiT(2*j+1)*p(2*Tdim+2);
                minusDS(2*Tdim+2) = minusDS(2*Tdim+2) + a*chiT(2*j+1)*p(2*j+1);
		    	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//boundaries			
			if (t==(NT-1))
				{
				kineticS += Dt*pow(Cp(neigh(j,1,1,Nb))-Cp(j),2.0)/a/2.0; //n.b. no contribution from time derivative term at the final time boundary
				pot_l += Dt*a*Va(Cp(j));
				pot_e += Dt*a*Vb(Cp(j));
				
				DDS.insert(2*j,2*(j-1)+1) = 1.0; //zero imaginary part of time  at final time boundary
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
			else if (t==0)
				{
				kineticS += Dt*pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
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
				kineticS += Dt*pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
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
                        unsigned int neighb = neigh(j,direc,sign,NT);
                        
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
	
        	break;
        } //ending while loop

} //ending parameter loop

return 0;
}

