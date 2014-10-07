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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//getting data from files

//getting mainInputs
unsigned long long int minFile, maxFile;
unsigned int loops;
double minTheta, maxTheta;
ifstream fmainin;
fmainin.open("mainInputs", ios::in);
string lastLine = getLastLine(fmainin);
istringstream ss(lastLine);
ss >> minFile >> maxFile >> minTheta >> maxTheta >> loops;
fmainin.close();

//defining the timeNumber
string timeNumber = currentDateTime();

//copying a version of mainInputs with timeNumber
string runInputs = "./data/" + timeNumber + "mainInputs";
copyFile("mainInputs",runInputs);

//getting list of relevant data files
system("dir ./data/* > dataFiles");
vector<string> filenames, piFiles, inputsFiles, eigenvectorFiles, eigenvalueFiles;
filenames = readDataFiles(minFile,maxFile);
piFiles = findStrings(filenames,"tpip");
inputsFiles = findStrings(filenames,"inputs");
eigenvectorFiles = findStrings(filenames,"eigVec");
eigenvalueFiles = findStrings(filenames,"eigValue");
if (piFiles.size()!=inputsFiles.size() || piFiles.size()!=eigenvectorFiles.size() || piFiles.size()!=eigenvalueFiles.size() || piFiles.size()==0)
	{
	cout << "required files not available" << endl;
	return 0;
	}
vector <unsigned long long int> fileNumbers;
fileNumbers = getInts(piFiles);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//beginning file loop
for (unsigned int fileLoop=0; fileLoop<piFiles.size(); fileLoop++)
	{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//getting variables from inputs

	ifstream fin;
	fin.open(inputsFiles[fileLoop]);
	string line;
	unsigned int lineNumber = 0;
	while(getline(fin,line))
		{
		if(line[0] == '#')
			{
			continue;
			}
		if (lineNumber==0)
			{
			istringstream ss(line);
			ss >> N >> Na >> Nb >> Nc >> dE >> Tb >> theta;
			if (absolute(theta-minTheta)>2.0e-16)
				{
				cout << "minTheta != theta" << endl;
				cout << minTheta << " != " << theta << endl;
				}
			lineNumber++;
			}
		else if (lineNumber==1)
			{
			istringstream ss(line);
			ss >> aq.inputChoice >> aq.fileNo >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
			ss >> alpha >> open >> amp;
			lineNumber++;
			}
		else if(lineNumber==2)
			{
			istringstream ss(line);
			double temp;
			ss >> temp >> temp >> negEigDone;
			lineNumber++;
			}
		}
	fin.close();
	inP = aq.inputChoice;

	//derived quantities
	NT = Na + Nb + Nc;
	epsilon = dE;
	R = 2.0/3.0/epsilon;
	alpha *= R;
	Gamma = exp(-theta);
	vec negVec(2*N*Nb+1);
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
		//cout << "Tb>R so using negEig, need to have run pi with inP='b'" << endl;
		if (negEigDone==0)
			{
			system("./negEig");
			cout << "negEig run" << endl;
			}
		negVec = loadVector(eigenvectorFiles[fileLoop],Nb,1);
		ifstream eigFile;
		eigFile.open(eigenvalueFiles[fileLoop]);
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
	closenessE = 1.0e-2;
	closenessL = 1.0e-2;
	closenessT = 1.0e-5;
	closenessP = 0.5;

	string loop_choice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
	string print_choice = aq.printChoice;
	
	//defining some important scalar quantities
	double S1 = 2.0/3.0; //mass of kink multiplied by lambda
	double twaction = -pi*epsilon*pow(R,2)/2.0 + pi*R*S1;
	
	//finding minima of potential. solving p^3 + 0*p^2 + b*p + c = 0
	double b_parameter = -1.0;
	double c_parameter = -epsilon;
	gsl_poly_solve_cubic (0, b_parameter, c_parameter, &root[0], &root[1], &root[2]);
	sort(root.begin(),root.end());
	comp ergZero = N*a*V(root[0]);
	mass2 = real(ddV(root[0]));

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
	//beginning theta loop
	for (unsigned int loop=0; loop<loops; loop++)
		{
		if (loops>1)
			{
			theta = minTheta + (maxTheta - minTheta)*loop/(loops-1.0);
			Gamma = exp(-theta);
			}
		
		//defining a time and starting the clock
		clock_t time;
		time = clock();
	
		//printing loop name and parameters
		printf("%12s%12s\n","timeNumber: ",timeNumber.c_str());		
		printParameters();
	
		//defining energy and number vectors
		cVec erg(NT);
		vec linErg(NT);
		vec linNum(NT);
		
		//defining the action and bound and W
		comp action = twaction;
		double bound = 0.0;
		double W = 0.0;
		double E;
		double Num;
		
		//defining some quantities used to stop the Newton-Raphson loop when action stops varying
		comp action_last = action;
		unsigned int runs_count = 0;
		unsigned int min_runs = 2;
		vector<double> action_test(1);	action_test[0] = 1.0;
		vector<double> sol_test(1);	sol_test[0] = 1.0;
		vector<double> solM_test(1);	solM_test[0] = 1.0;
		vector<double> delta_test(1); delta_test[0] = 1.0;
		vector<double> calc_test(1); calc_test[0] = 1.0;
		vector<double> erg_test(1); erg_test[0] = 1.0;
		vector<double> lin_test(1); lin_test[0] = 1.0;
		vector<double> true_test(1); true_test[0] = 1.0;
		vector<double> mom_test(1); mom_test[0] = 1.0;

		//initializing phi (=p)
		vec p(2*N*NT+2);
		if (loop==0)
			{
			p = loadVector(piFiles[fileLoop],NT,2);
			}
		else
			{
			unsigned int toLoad = loop-1;
			string loadfile = "./data/" + timeNumber + "main" + to_string(toLoad)+".dat";
			p = loadVector(loadfile,NT,2);
			}
			
		//very early vector print
		string earlyPrintFile = "data/" + timeNumber + "mainE"+ to_string(fileLoop) + to_string(loop) + "0.dat";
		printVectorB(earlyPrintFile,p);
		
		//defining complexified vector Cp
		cVec Cp(NT*N);
		Cp = vecComplex(p,N*NT);
	
		//defining DDS and minusDS
		spMat DDS(2*N*NT+2,2*N*NT+2);
		vec minusDS(2*N*NT+2);
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//beginning newton-raphson loop	
		while ((sol_test.back()>closenessS || solM_test.back()>closenessSM || runs_count<min_runs))
			{
			runs_count++;
			
			//defining the zero mode at the final time boundary and the time step before
			vec chiX(NT*N);	chiX = Eigen::VectorXd::Zero(N*NT); //to fix spatial zero mode
			vec chiT(NT*N);	chiT = Eigen::VectorXd::Zero(N*NT); //to fix real time zero mode
			for (unsigned int j=0; j<N; j++)
				{
				unsigned int posC = j*NT+(Na+Nb-1);
				unsigned int posA = j*NT;
				unsigned int posDm = j*NT+(NT-2);
				chiX(posC) = p(2*neigh(posC,1,1,NT))-p(2*neigh(posC,1,-1,NT));
				chiT(posA) = p(2*(posA+1))-p(2*posA); //could also fix against negVec if this doesn't work
				chiT(posDm) = p(2*(posDm+1))-p(2*posDm);
				}
			
			// allocating memory for DS, DDS
			minusDS = Eigen::VectorXd::Zero(2*N*NT+2); //initializing to zero
			DDS.setZero(); //just making sure
			Eigen::VectorXi DDS_to_reserve(2*N*NT+2);//number of non-zero elements per column
			DDS_to_reserve = Eigen::VectorXi::Constant(2*N*NT+2,11);
			DDS_to_reserve(0) = 3; //these need to be changed when boundary conditions need to be more compicated
			DDS_to_reserve(1) = 3;
			DDS_to_reserve(2*N*NT-2) = 3;
			DDS_to_reserve(2*N*NT-1) = 3;
			DDS_to_reserve(2*N*NT) = N;
			DDS_to_reserve(2*N*NT+1) = 2*N;
			DDS.reserve(DDS_to_reserve);
			
			//initializing to zero
			comp kineticS = 0.0;
			comp kineticT = 0.0;
			comp pot_l = 0.0;
			comp pot_e = 0.0;
			bound = 0.0;
			erg = Eigen::VectorXcd::Constant(NT,-ergZero);
			linErg = Eigen::VectorXd::Zero(NT);
			linNum = Eigen::VectorXd::Zero(NT);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//assigning values to minusDS and DDS and evaluating action
			for (unsigned long int j = 0; j < N*NT; j++)
				{		
				unsigned int t = intCoord(j,0,NT); //coordinates
				unsigned int x = intCoord(j,1,NT);
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
					DDS.insert(2*(j+1),2*N*NT+1) = a*chiT(j); //there should be no chiT on the final time slice or this line will go wrong
					DDS.insert(2*N*NT+1,2*(j+1)) = a*chiT(j);
					DDS.insert(2*j,2*N*NT+1) = -a*chiT(j);
					DDS.insert(2*N*NT+1,2*j) = -a*chiT(j);
		            minusDS(2*(j+1)) += - a*chiT(j)*p(2*N*NT+1);
		            minusDS(2*N*NT+1) += - a*chiT(j)*p(2*j);
		            minusDS(2*j) += a*chiT(j)*p(2*N*NT+1);
		            minusDS(2*N*NT+1) += a*chiT(j)*p(2*j);
					}
					
				if (absolute(theta)<2.0e-16)
					{
					for (unsigned int k=0;k<N;k++)
						{
						unsigned int l = k*NT+t;
						linErg(t) += Eomega(x,k)*(p(2*l)-root[0])*(p(2*j)-root[0]) + Eomega(x,k)*p(2*j+1)*p(2*l+1); //middle sign may be negative - check this
						linNum(t) += omega(x,k)*(p(2*l)-root[0])*(p(2*j)-root[0]) + omega(x,k)*p(2*j+1)*p(2*l+1);
						}
					}
				else
					{
					for (unsigned int k=0;k<N;k++)
						{
						unsigned int l = k*NT+t;
						linErg(t) += 2.0*Gamma*omega(x,k)*(p(2*l)-root[0])*(p(2*j)-root[0])/pow(1.0+Gamma,2.0) + 2.0*Gamma*omega(x,k)*p(2*j+1)*p(2*l+1)/pow(1.0-Gamma,2.0);
						linNum(t) += 2.0*Gamma*Eomega(x,k)*(p(2*l)-root[0])*(p(2*j)-root[0])/pow(1.0+Gamma,2.0) + 2.0*Gamma*Eomega(x,k)*p(2*j+1)*p(2*l+1)/pow(1.0-Gamma,2.0);
						}
					}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//boundaries			
				if (t==(NT-1))
					{
					comp tempErg =  (kineticS + pot_l + pot_e)/Dt;
					kineticS += Dt*pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0;
					pot_l += Dt*a*Va(Cp(j));
					pot_e += Dt*a*Vb(Cp(j));
					erg(t) += + (kineticS + pot_l + pot_e)/Dt - tempErg;
				
					DDS.insert(2*j,2*(j-1)+1) = 1.0; //zero imaginary part of time derivative
					DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
					}
				else if (t==0)
					{
					comp tempErg =  (kineticS + pot_l + pot_e)/Dt + kineticT/dt;
					kineticS += Dt*pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0;
					kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					pot_l += Dt*a*Va(Cp(j));
					pot_e += Dt*a*Vb(Cp(j));
					erg(t) += (kineticS + pot_l + pot_e)/Dt + kineticT/dt - tempErg;
					
					if (absolute(theta)<2.0e-16)
						{
						DDS.insert(2*j,2*(j+1)+1) = 1.0; //zero imaginary part of time derivative
						DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
						}
					else
						{
						for (unsigned int k=1; k<2*2; k++)
				        	{
				            int sign = pow(-1,k+1);
				            int direc = (int)(k/2.0);
				            if (direc == 0)
				            	{
				                minusDS(2*j) += real(a*Cp(j+sign)/dt);
				                minusDS(2*j+1) += imag(a*Cp(j+sign)/dt);
				                }
				            else
				            	{
				                unsigned int neighb = neigh(j,direc,sign,NT);
				                minusDS(2*j) += - real(Dt*Cp(neighb)/a);
				                minusDS(2*j+1) += - imag(Dt*Cp(neighb)/a);
				                }
				            }
                        comp temp0 = a/dt;
		            	comp temp1 = a*Dt*(2.0*Cp(j)/pow(a,2.0));//dropped this last term, as linearised + dV(Cp(j)));

				        minusDS(2*j) += real(temp1 - temp0*Cp(j));
				        minusDS(2*j+1) += imag(temp1 - temp0*Cp(j));

                        for (unsigned int k=0; k<N; k++)
                        	{
                        	//dropped two terms here, as had imag(real()) or vice versa
                        	DDS.insert(2*j,2*k*NT+1) = -omega(x,k)*(1.0-Gamma)/(1.0+Gamma);
                        	DDS.insert(2*j+1,2*k*NT) = omega(x,k)*(1.0+Gamma)/(1.0-Gamma);
							unsigned int m = x*NT;
							unsigned int n = k*NT;
							bound += -(1.0-Gamma)*omega(m,n)*(p(2*m)-root[0])*(p(2*n)-root[0])/(1.0+Gamma) + (1.0+Gamma)*omega(m,n)*p(2*m+1)*p(2*n+1)/(1.0-Gamma);
                        	}
						}
					}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//bulk
				else
					{
					comp tempErg =  (kineticS + pot_l + pot_e)/Dt + kineticT/dt;
					kineticS += Dt*pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0;
					kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					pot_l += Dt*a*Va(Cp(j));
					pot_e += Dt*a*Vb(Cp(j));
					erg(t) += (kineticS + pot_l + pot_e)/Dt + kineticT/dt - tempErg;
				
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
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			if ((print_choice.compare("a")==0 || print_choice.compare("e")==0) && 1==0) //have stopped this one as it's annoying
				{
				printAction(kineticT-kineticS,pot_l,pot_e);
				}
			if ((print_choice.compare("v")==0 || print_choice.compare("e")==0) && 1==0)
				{
				string minusDSfile = "./data/" + timeNumber + "mainminusDSE"+to_string(fileLoop)+to_string(loop)+to_string(runs_count)+".dat";
				printVector(minusDSfile,minusDS);
				}
			if ((print_choice.compare("p")==0 || print_choice.compare("e")==0) && delta_test.back()>0.2)
				{
				string piEarlyFile = "./data/" + timeNumber + "mainpE"+to_string(fileLoop)+to_string(loop)+to_string(runs_count)+".dat";
				printVector(piEarlyFile,p);
				}
			if ((print_choice.compare("m")==0 || print_choice.compare("e")==0) && 1==0)
				{
				string DDSfile = "./data/" + timeNumber + "mainDDSE"+to_string(fileLoop)+to_string(loop)+to_string(runs_count)+".dat";
				printSpmat(DDSfile,DDS);
				}
			}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			//checking linErg, linNum, bound, computing W, E
			
			//checking linearisation of linErg and linNum
			double linTestE;	double linTestN;
			double linEMax = 0.0;	double linEMin = 1.0e6; //surely it's going to be less than this
			double linNMax = 0.0;	double linNMin = 1.0e6;
			unsigned int linearInt = (int)(Na/6);
			for (unsigned int j=1;j<(linearInt+1);j++)
				{
				if (absolute(linErg(j))>linEMax)
					{
					linEMax = absolute(linErg(j));
					}
				if (absolute(linErg(j))<linEMin)
					{
					linEMin = absolute(linErg(j));
					}
				if (absolute(linNum(j))>linNMax)
					{
					linNMax = absolute(linNum(j));
					}
				if (absolute(linNum(j))<linNMin)
					{
					linNMin = absolute(linNum(j));
					}
				}
			linTestE = (linEMax-linEMin)*2.0/(linEMax+linEMin);
			linTestN = (linNMax-linNMin)*2.0/(linNMax+linNMin);
			lin_test.push_back(linTestE);
			if (linTestN>linTestE)
				{
				lin_test.back() = linTestN;
				}
				
			//checking agreement between erg and linErg
			double trueTest = real(erg(1))-linErg(1); //not using zero as boundaries are funny
			trueTest = trueTest*2.0/(real(erg(1))+linErg(1));
			trueTest = absolute(trueTest);
			true_test.push_back(trueTest);
			
			//checking conservation of E
			double ergTest = real(erg(1)-erg(NT-2));
			ergTest = ergTest*2.0/(real(erg(1)+erg(NT-2)));
			ergTest = absolute(ergTest);
			erg_test.push_back(ergTest);
			if (ergTest>closenessE)
				{
				cout << "ergTest = " << ergTest << endl;
				}
							
			//checking lattice small enough for E, should have parameter for this
			double momTest = E*a/Num*pi;
			mom_test.push_back(momTest);
			if (momTest>closenessP)
				{
				cout << "lattice not small enough for momenta" << endl;
				cout << "momTest = " << momTest << endl;
				}
						
			//defining E, Num and cW
			E = real(erg(1));
			Num = linNum(1);
			W = - E*2.0*Tb - theta*Num - bound + 2.0*imag(action);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

			//solving for delta in DDS*delta=minusDS, where p' = p + delta		
			vec delta(2*N*NT+2);
			delta = Eigen::VectorXd::Zero(2*N*NT+2);
			DDS.makeCompressed();
			Eigen::SparseLU<spMat> solver;
			break;
		
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
			vec diff(2*N*NT+2);
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
			Cp = vecComplex(p,N*NT);
		
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
				printf("%6s%6s%14s%14s%14s%14s%14s%14s\n","loop","runsCount","actionTest","solTest","solMTest","deltaTest","linTest","trueTest");
				}
			printf("%6i%6i%14g%14g%14g%14g%14g%14g\n",loop,runs_count,action_test.back(),sol_test.back(),solM_test.back(),delta_test.back(),lin_test.back(),true_test.back());
			
			break;
			} //ending while loop
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	    		//misc end of program tasks - mostly printing

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
		
		//stopping clock
		time = clock() - time;
		double realtime = time/1000000.0;
	
		//printing to terminal
		printf("\n");
		printf("%8s%8s%8s%8s%8s%8s%8s%8s%16s%16s%16s%16s\n","runs","time","N","NT","L","Tb","dE","theta","E","im(action)","bound","W");
		printf("%8i%8g%8i%8i%8g%8g%8g%8g%16g%16g%16g%16g\n",runs_count,realtime,N,NT,L,Tb,dE,theta,E,imag(action),bound,real(W));
		printf("\n");
		 printf("%60s\n","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

		//printing action value
		FILE * actionfile;
		actionfile = fopen("./data/mainAction.dat","a");
		fprintf(actionfile,"%16s%8i%8i%8g%8g%8g%8g%14g%12g%14g%14g%14g%14g\n",timeNumber.c_str(),N,NT,L,Tb,dE,theta,E\
		,imag(action),real(W),sol_test.back(),solM_test.back(),lin_test.back());
		fclose(actionfile);
	
		//copying a version of inputs with timeNumber and theta changed
		string runInputs = "./data/" + timeNumber + "inputs";
		changeInputs(runInputs, "theta", to_string(theta));
	
		//printing output phi on whole time contour
		string tpifile = "./data/" + timeNumber + "mainp"+to_string(fileLoop)+to_string(loop)+".dat";
		printVector(tpifile,p);
		gp(tpifile,"repi.gp");
	
		//printing output minusDS				
		string minusDSfile = "./data/" + timeNumber + "mainminusDS"+to_string(fileLoop)+to_string(loop)+".dat";
		printVector(minusDSfile,minusDS);
				
		//printing output DDS
		string DDSfile = "./data/" + timeNumber + "mainDDS"+to_string(fileLoop)+to_string(loop)+".dat";
		printSpmat(DDSfile,DDS);
	
		//printing linErgVec
		//string linErgFile = "./data/" + timeNumber + "mainlinErg"+to_string(fileLoop)+to_string(loop)+".dat";
		//simplePrintCVector(linErgFile,linErgA);
		//gpSimple(linErgFile);
	
		//printing erg
		string ergFile = "./data/" + timeNumber + "erg" + to_string(fileLoop)+to_string(loop)+".dat";
		simplePrintCVector(ergFile,erg);
		//gpSimple(ergFile);
		
		} //ending theta loop
	} //ending file loop

return 0;
}

