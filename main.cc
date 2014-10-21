#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cmath>
#include "pf.h"
#include "files.h"

int main()
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//getting data from files

//getting mainInputs
unsigned long long int minFile, maxFile;
unsigned int loops;
double minTb, maxTb, minTheta, maxTheta;
ifstream fmainin;
fmainin.open("mainInputs", ios::in);
if (fmainin.is_open())
	{
	string line;
	unsigned int lineNumber = 0;
	while(getline(fmainin,line))
		{
		if(line[0] == '#')
			{
			continue;
			}
		istringstream ss(line);
		if (lineNumber==0)
			{
			ss >> inF >> minFile >> maxFile >> firstLoop;
			lineNumber++;
			}
		else if (lineNumber==1)
			{
			ss >> minTb >> maxTb >> minTheta >> maxTheta >> loops >> zmt >> zmx;
			lineNumber++;
			}
		}
	}
else
	{
	cout << "unable to open mainInputs" << endl;
	}
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
if (inF.compare("p")==0)
	{
	piFiles = findStrings(filenames,"tpip");
	inputsFiles = findStrings(filenames,"inputsPi");
	eigenvectorFiles = findStrings(filenames,"eigVec");
	eigenvalueFiles = findStrings(filenames,"eigValue");
	if  (eigenvectorFiles.size()!=1 || eigenvalueFiles.size()!=1)
		{
		cout << "correct eigen files not available" << endl;
		cout << "eigenvalueFiles.size() = " << eigenvalueFiles.size() << endl;
		cout << "eigenvectorFiles.size() = " << eigenvectorFiles.size() << endl;
		}
	}
else if (inF.compare("m")==0)
	{
	piFiles = findStrings(filenames,"mainp");
	inputsFiles = findStrings(filenames,"inputsM");
	}
else
	{
	cout << "inF error" << endl;
	}
vector <vector<string>*> files = {&piFiles,&inputsFiles};
for (unsigned int k=0;k<files.size();k++)
	{
	vector <string>* tempVecStr = files[k];
	vector <unsigned int> loopNumbers = getLastInts(*tempVecStr);
	for (unsigned int l=0;l<(*tempVecStr).size();l++)
		{
		if (loopNumbers[l]<firstLoop)
			{
			(*tempVecStr).erase((*tempVecStr).begin()+l);
			}
		}
	}
if (piFiles.size()!=inputsFiles.size() || piFiles.size()==0)
	{
	for (unsigned int j=0;j<files.size();j++)
		{
		for (unsigned int k=0;k<files.size();k++)
			{
			*files[j] = reduceTo(*files[j],*files[k]);
			}
		}
	if (piFiles.size()!=inputsFiles.size() || piFiles.size()==0)
		{
		cout << "required files not available" << endl;
		cout << "piFiles.size() = " << piFiles.size() << endl;
		cout << "inputsFiles.size() = " << inputsFiles.size() << endl;
		return 0;
		}
	}
		
sort(piFiles.begin(), piFiles.end());
sort(inputsFiles.begin(), inputsFiles.end());

vector <unsigned long long int> fileNumbers = getInts(piFiles);
cout << endl;
for (unsigned int j=0; j<inputsFiles.size();j++)
	{
	cout << inputsFiles[j] << " " << piFiles[j] << endl;
	}
cout << endl;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//beginning file loop
for (unsigned int fileLoop=0; fileLoop<piFiles.size(); fileLoop++)
	{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//getting variables from inputs

	ifstream fin;
	fin.open(inputsFiles[fileLoop]);
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
				if (absolute(theta-minTheta)>2.0e-16 && loops>1)
					{
					cout << "minTheta != theta" << endl;
					cout << minTheta << " != " << theta << endl;
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
		cout << "unable to open " << inputsFiles[fileLoop] << endl;
		}
	fin.close();
	inP = aq.inputChoice;
	
	//potential functions
	if (pot[0]=='1')
		{
		V = &V1;
		V0 = &V10;
		Ve = &V1e;
		dV = &dV1;
		ddV = &ddV1;
		}
	else if (pot[0]=='2')
		{
		V = &V2;
		V0 = &V20;
		Ve = &V2e;
		dV = &dV2;
		ddV = &ddV2;
		}

	//derived quantities
	NT = Na + Nb + Nc;
	epsilon = dE;
	R = 2.0/3.0/epsilon;
	//alpha *= R;
	Gamma = exp(-theta);
	vec negVec(2*N*Nb+1);
	L = LoR*R;
	if (Tb<R)
		{
		angle = asin(Tb/R);
		double Ltemp = 1.5*(1.5*Tb*tan(angle));
		if (Ltemp<L) //making sure to use the smaller of the two possible Ls
			{
			L=Ltemp;
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
	closenessC = 1.0e-12;
	closenessE = 1.0e-2;
	closenessL = 1.0e-2;
	closenessT = 1.0e-5;
	closenessP = 0.5;
	closenessR = 1.0e-4;

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

	if (zmt[0]=='n' || zmx[0]=='n')
		{
		//cout << "Tb>R so using negEig, need to have run pi with inP='b'" << endl;
		if (negEigDone==0)
			{
			system("./negEig"); //negEig now needs timeNumber
			char * fileNumber = (char *)(fileNumbers[fileLoop]);
			system(fileNumber); //won't work
			cout << "negEig run" << endl;
			}
		negVec = loadVector(eigenvectorFiles[0],Nb,1);
		ifstream eigFile;
		eigFile.open(eigenvalueFiles[0]);
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//beginning theta loop
	for (unsigned int loop=0; loop<loops; loop++)
		{
		if (loops>1)
			{
			if (absolute(theta-minTheta)>2.0e-16 && absolute(Tb-minTb)>2.0e-16 && loop==0)
				{
				cout << "inputs and mainInputs don't agree" << endl;
				}
			if (absolute(maxTheta-minTheta)>2.0e-16)
				{
				theta = minTheta + (maxTheta - minTheta)*loop/(loops-1.0);
				Gamma = exp(-theta);
				}
			else if (absolute(maxTb-minTb)>2.0e-16)
				{
				Tb = minTb + (maxTb - minTb)*loop/(loops-1.0);
				changeDouble ("Tb",Tb);
				}
			else
				{
				cout << "nothing to loop over, set loops to 1" << endl;
				}
			}
		
		//defining a time and starting the clock
		clock_t time;
		time = clock();
	
		//defining energy and number vectors
		cVec erg(NT);
		vec linErg(NT);
		vec linNum(NT);
		
		//defining the action and bound and W
		comp action = i*twaction;
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
		vector<double> reg_test(1); reg_test[0] = 1.0;
		
		//printing loop name and parameters
		printf("%12s%12s\n","timeNumber: ",timeNumber.c_str());

		//initializing phi (=p)
		vec p(2*N*NT+2);
		if (loop==0)
			{
			p = loadVector(piFiles[fileLoop],NT,2);
			printf("%12s%12llu\n","input: ",fileNumbers[fileLoop]);
			}
		else
			{
			unsigned int toLoad = loop-1;
			string loadfile = "./data/" + timeNumber + "mainp" + to_string(fileLoop) + to_string(toLoad)+".dat";
			p = loadVector(loadfile,NT,2);
			printf("%12s%12s%12s%12u\n","input: ",timeNumber.c_str(), ", loop: ", toLoad);
			}
				
		printParameters();
			
		//very early vector print
		string earlyPrintFile = "data/" + timeNumber + "mainE"+ to_string(fileLoop) + "_" + to_string(loop) + "_" + "0.dat";
		printVector(earlyPrintFile,p);
		
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
				unsigned int posX, posT, posCe;
				char* tempChar = &zmx.back();
				stringstream ssX(tempChar);
				unsigned int slicesX, slicesT;
				ssX >> slicesX;
				tempChar = &zmt.back();
				stringstream ssT(tempChar);
				ssT >> slicesT; //all seems v long winded but other ways seemed to fail
				posCe = j*Nb+Nb-slicesT; //position C for Euclidean vector, i.e. for negVec
				map<char,unsigned int> posMap;
				posMap['A'] = j*NT;
				posMap['B'] = j*NT + (Na-1);
				posMap['C'] = j*NT+Na+Nb-slicesT;
				posMap['D'] = j*NT+(NT-1)-slicesT;
				if (zmt.size()<3) { cout << "zmt lacks info, zmt = " << zmt << endl; }
				for (unsigned int l=0;l<(zmt.size()-2);l++)
					{
					if (posMap.find(zmt[1+l])!=posMap.end())
						{
						posT = posMap.at(zmt[1+l]);
						for (unsigned int k=0;k<slicesT;k++)
							{
							if (zmt[0]=='n')
								{
								chiT(posT+k) = negVec(2*(posCe+k));
								}
							else if (zmt[0]=='d')
								{
								chiT(posT+k) = p(2*(posT+k+1))-p(2*(posT+k));
								}
							else
								{
								cout << "choice of zmt not allowed" << endl;
								}
							}
						}
					}
				posMap.erase('C');
				posMap.erase('D');
				posMap['C'] = j*NT+Na+Nb-slicesX;
				posMap['D'] = j*NT+NT-slicesX;
				posCe = j*Nb+Nb-slicesX;
				if (zmx.size()<3) { cout << "zmx lacks info, zmx = " << zmx << endl; }
				for (unsigned int l=0;l<(zmx.size()-2);l++)
					{
					if (posMap.find(zmx[1+l])!=posMap.end())
						{
						posX = posMap.at(zmx[1+l]);
						for (unsigned int k=0;k<slicesX;k++)
							{
							if (zmx[0]=='n')
								{
								chiT(posX+k) = negVec(2*(posCe+k));
								}
							else if (zmx[0]=='d')
								{
								chiX(posX+k) = p(2*neigh(posX+k,1,1,NT))-p(2*neigh(posX+k,1,-1,NT));
								}
							else
								{
								cout << "choice of zmx not allowed" << endl;
								}
							}
						}
					}
				}
			double normX = chiX.dot(chiX);
			normX = pow(normX,0.5);
			double normT = chiT.dot(chiT);
			normT = pow(normT,0.5);
			if (absolute(normX)<2.0e-16 || absolute(normT)<2.0e-16)
				{
				cout << "norm of chiX = " << normX << ", norm of chiT = " << normT << endl;
				}
			chiX = chiX/normX;
			chiT = chiT/normT;
			
			// allocating memory for DS, DDS
			minusDS = Eigen::VectorXd::Zero(2*N*NT+2); //initializing to zero
			DDS.setZero(); //just making sure
			Eigen::VectorXi DDS_to_reserve(2*N*NT+2);//number of non-zero elements per column
			DDS_to_reserve = Eigen::VectorXi::Constant(2*N*NT+2,11);
			if (absolute(theta)<2.0e-16)
				{
				DDS_to_reserve(0) = 3;
				DDS_to_reserve(1) = 3;
				}
			else
				{
				DDS_to_reserve(0) = N+1;
				DDS_to_reserve(1) = N+1;
				}
			DDS_to_reserve(2*N*NT-2) = 3;
			DDS_to_reserve(2*N*NT-1) = 3;
			DDS_to_reserve(2*N*NT) = N*(zmx.size()-2);
			DDS_to_reserve(2*N*NT+1) = 2*N*(zmt.size()-2);
			DDS.reserve(DDS_to_reserve);
			
			//initializing to zero
			comp kineticS = 0.0;
			comp kineticT = 0.0;
			comp pot_l = 0.0;
			comp pot_e = 0.0;
			comp pot_r = 0.0;
			bound = 0.0;
			erg = Eigen::VectorXcd::Constant(NT,-ergZero);
			linErg = Eigen::VectorXd::Zero(NT);
			linNum = Eigen::VectorXd::Zero(NT);
			
			//lambda functions for pot_r
			auto Vr = [&] (const comp & phi)
				{
				return reg*VrFn(phi,root[0],root[2]);
				};
			auto dVr = [&] (const comp & phi)
				{
				return reg*dVrFn(phi,root[0],root[2]);
				};
			auto ddVr = [&] (const comp & phi)
				{
				return reg*ddVrFn(phi,root[0],root[2]);
				};

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
					DDS.coeffRef(2*(j+1),2*N*NT+1) += a*chiT(j); //chiT should be 0 at t=(NT-1) or this line will go wrong
					DDS.coeffRef(2*N*NT+1,2*(j+1)) += a*chiT(j);
					DDS.coeffRef(2*j,2*N*NT+1) += -a*chiT(j);
					DDS.coeffRef(2*N*NT+1,2*j) += -a*chiT(j);
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
						linNum(t) += 2.0*Gamma*omega(x,k)*(p(2*l)-root[0])*(p(2*j)-root[0])/pow(1.0+Gamma,2.0) - 2.0*Gamma*omega(x,k)*p(2*j+1)*p(2*l+1)/pow(1.0-Gamma,2.0); //are signs right? shouldn't this be positive definite?
						linErg(t) += 2.0*Gamma*Eomega(x,k)*(p(2*l)-root[0])*(p(2*j)-root[0])/pow(1.0+Gamma,2.0) - 2.0*Gamma*Eomega(x,k)*p(2*j+1)*p(2*l+1)/pow(1.0-Gamma,2.0);
						}
					}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//boundaries			
				if (t==(NT-1))
					{
					kineticS += Dt*pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0;
					pot_l += Dt*a*V0(Cp(j));
					pot_e += Dt*a*Ve(Cp(j));
					pot_r += Dt*a*Vr(Cp(j));
					erg(t) += pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
					DDS.insert(2*j,2*(j-1)+1) = 1.0; //zero imaginary part of time derivative
					DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
					}
				else if (t==0)
					{
					kineticS += Dt*pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0;
					kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					pot_l += Dt*a*V0(Cp(j));
					pot_e += Dt*a*Ve(Cp(j));
					pot_r += Dt*a*Vr(Cp(j));
					erg(t) += a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
					
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
		            	comp temp1 = a*Dt*(2.0*Cp(j)/pow(a,2.0));//dropped this last term, as linearised + dV(Cp(j))) + dVr(Cp(j));

				        minusDS(2*j) += real(temp1 - temp0*Cp(j));
				        minusDS(2*j+1) += imag(temp1 - temp0*Cp(j));

                        for (unsigned int k=0; k<N; k++)
                        	{
                        	//dropped two terms here, as had imag(real()) or vice versa
							unsigned int m = x*NT;
							unsigned int n = k*NT;
							DDS.insert(2*j,2*n+1) = -omega(x,k)*(1.0-Gamma)/(1.0+Gamma);
                        	DDS.insert(2*j+1,2*n) = omega(x,k)*(1.0+Gamma)/(1.0-Gamma);
							bound += -(1.0-Gamma)*omega(x,k)*(p(2*m)-root[0])*(p(2*n)-root[0])/(1.0+Gamma) + (1.0+Gamma)*omega(x,k)*p(2*m+1)*p(2*n+1)/(1.0-Gamma);
                        	}
						}
					}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//bulk
				else
					{
					kineticS += Dt*pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0;
					kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					pot_l += Dt*a*V0(Cp(j));
					pot_e += Dt*a*Ve(Cp(j));
					pot_r += Dt*a*Vr(Cp(j));
					erg(t) += a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neigh(j,1,1,NT))-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
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
			//checking linErg, linNum, bound, computing W, E
			
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
			
			//checking linearisation of linErg and linNum
			double linTestE;	double linTestN;
			double linEMax = 0.0;	double linEMin = 5.0e15; //surely it's going to be less than this
			double linNMax = 0.0;	double linNMin = 5.0e15;
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
						
			//defining E, Num and cW
			E = 0;
			Num = 0;
			for (unsigned int j=0; j<linearInt; j++)
				{
				E += real(linErg(j));
				Num += real(linNum(j));
				}
			E /= linearInt;
			Num /= linearInt;
			W = - E*2.0*Tb - theta*Num - bound + 2.0*imag(action);
			
			//checking lattice small enough for E, should have parameter for this
			double momTest = E*a/Num*pi;
			mom_test.push_back(momTest);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			string prefix = "./data/" + timeNumber;
			string suffix = to_string(fileLoop)+"_"+to_string(loop)+"_"+to_string(runs_count)+".dat";
			if ((print_choice.compare("a")==0 || print_choice.compare("e")==0) && 1==0)
				{
				printAction(kineticT-kineticS,pot_l,pot_e);
				}
			if ((print_choice.compare("v")==0 || print_choice.compare("e")==0))
				{
				string minusDSfile = prefix + "mainminusDSE"+suffix;
				printVector(minusDSfile,minusDS);
				}
			if ((print_choice.compare("p")==0 || print_choice.compare("e")==0) || delta_test.back()>0.2)
				{
				string piEarlyFile = prefix + "mainE"+suffix;
				printVector(piEarlyFile,p);
				//gp(piEarlyFile,"repi.gp");
				}
			if ((print_choice.compare("m")==0 || print_choice.compare("e")==0))
				{
				string DDSfile = prefix + "mainDDSE"+suffix;
				printSpmat(DDSfile,DDS);
				}
			if ((print_choice.compare("z")==0 || print_choice.compare("e")==0))
				{
				string earlychiXFile = prefix + "mainchiXE" + suffix;
				printVector(earlychiXFile,chiX);
				//gp(earlychiXFile,"repi.gp");
				string earlychiTFile = prefix + "mainchiTE" + suffix;
				printVector(earlychiTFile,chiT);
				//gp(earlychiTFile,"repi.gp");
				}
			if ((print_choice.compare("l")==0 || print_choice.compare("e")==0))
				{
				string earlyLinErgFile = prefix + "mainlinErgE"+suffix;
				simplePrintVector(earlyLinErgFile,linErg);
				//gpSimple(earlyLinErgFile);
				string earlyErgFile = prefix + "mainergE" + suffix;
				simplePrintCVector(earlyErgFile,erg);
				//gpSimple(earlyErgFile);
				}
			}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			//solving for delta in DDS*delta=minusDS, where p' = p + delta		
			vec delta(2*N*NT+2);
			delta = Eigen::VectorXd::Zero(2*N*NT+2);
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
				printf("%10s%10s%14s%14s%14s%14s%14s%14s\n","loop","runsCount","actionTest","solTest","solMTest","deltaTest","linTest","trueTest");
				}
			printf("%10i%10i%14g%14g%14g%14g%14g%14g\n",loop,runs_count,action_test.back(),sol_test.back(),solM_test.back(),delta_test.back(),lin_test.back(),true_test.back());
			
			} //ending while loop
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	    		//misc end of program tasks - mostly printing

		//checking energy conserved
		if (erg_test.back()>closenessE)
				{
				cout << "test of energy conservation, ergTest = " << erg_test.back() << endl;
				}
		
		//checking lattice small enough
		if (mom_test.back()>closenessP)
			{
			cout << "lattice not small enough for momenta, momTest = "<< mom_test.back()  << endl;
			}
		
		//stopping clock
		time = clock() - time;
		double realtime = time/1000000.0;
	
		//printing to terminal
		printf("\n");
		printf("%8s%8s%8s%8s%8s%8s%8s%8s%14s%14s%14s\n","runs","time","N","NT","L","Tb","dE","theta","Num","E","W");
		printf("%8i%8g%8i%8i%8g%8g%8g%8g%14g%14g%14g\n",runs_count,realtime,N,NT,L,Tb,dE,theta,Num,E,real(W));
		printf("\n");
		 printf("%60s\n","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

		//printing action value
		FILE * actionfile;
		actionfile = fopen("./data/mainAction.dat","a");
		fprintf(actionfile,"%14s%6i%6i%8g%8g%8g%6g%14g%14g%12g%14g%14g\n",timeNumber.c_str(),N,NT,L,Tb,dE,theta,E,Num\
		,real(W),sol_test.back(),lin_test.back());
		fclose(actionfile);
	
		//copying a version of inputs with timeNumber and theta changed
		string runInputs = "./data/" + timeNumber + "inputsM"+ to_string(fileLoop) + "_" + to_string(loop);
		if (absolute(maxTheta-minTheta)>2.0e-16)
			{
			changeInputs(runInputs, "theta", to_string(theta));
			}
		else if (absolute(maxTb-minTb)>2.0e-16)
			{
			changeInputs(runInputs, "Tb", to_string(Tb));
			}
		else
			{
			copyFile("inputs",runInputs);
			}
		
		string prefix = "./data/" + timeNumber;
		string suffix = to_string(fileLoop)+"_" + to_string(loop)+".dat";
	
		//printing output phi
		string tpifile =  prefix + "mainp"+suffix;
		printVector(tpifile,p);
		gp(tpifile,"repi.gp");
	
		//printing output minusDS				
		string minusDSfile = prefix + "mainminusDS"+suffix;
		printVector(minusDSfile,minusDS);
				
		//printing output DDS
		string DDSfile = prefix + "mainDDS"+suffix;
		printSpmat(DDSfile,DDS);
	
		//printing linNumVec
		string linNumFile = prefix + "mainlinNum"+suffix;
		linNum.conservativeResize(Na);
		simplePrintVector(linNumFile,linNum);
		//gpSimple(linNumFile);	
	
		//printing linErgVec
		string linErgFile = prefix + "mainlinErg"+suffix;
		linErg.conservativeResize(Na);
		simplePrintVector(linErgFile,linErg);
		//gpSimple(linErgFile);
	
		//printing erg
		string ergFile = prefix + "mainerg" + suffix;
		erg.conservativeResize(Na);
		simplePrintCVector(ergFile,erg);
		//gpSimple(ergFile);
		
		} //ending theta loop
	} //ending file loop

return 0;
}

