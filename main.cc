#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cmath>
#include <limits>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include "pf.h"
#include "files.h"

int main()
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//getting data from mainInputs

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//sorting out files to load

//getting list of relevant data files, with timeNumbers between minFile and maxFile
int systemCall = system("dir ./data/* > dataFiles");
if (systemCall==-1)
	{
	cout << "system call failure, finding dataFiles" << endl;
	}
vector<string> filenames, piFiles, inputsFiles, eigenvectorFiles, eigenvalueFiles;
filenames = readDataFiles(minFile,maxFile);

//picking subset based on core of filename, using inF
if (inF.compare("p")==0)
	{
	piFiles = findStrings(filenames,"tpip");
	inputsFiles = findStrings(filenames,"inputsPi");
	eigenvectorFiles = findStrings(filenames,"eigVec");
	eigenvalueFiles = findStrings(filenames,"eigValue");
	if  (eigenvectorFiles.size()==0 || eigenvalueFiles.size()==0)
		{
		cout << "eigen files not available" << endl;
		cout << "eigenvalueFiles.size() = " << eigenvalueFiles.size() << endl;
		cout << "eigenvectorFiles.size() = " << eigenvectorFiles.size() << endl;
		}
	}
else if (inF.compare("m")==0)
	{
	piFiles = findStrings(filenames,"mainpi_");
	inputsFiles = findStrings(filenames,"inputsM");
	}
else
	{
	cout << "inF error" << endl;
	}

//picking subset based on loopNumbers being above firstLoop
vector <vector<string>*> files;
files.reserve(2);
files.push_back(&piFiles); files.push_back(&inputsFiles);
if (files.size()==0)
	{
	cout << "no files found" << endl;
	}
for (unsigned int k=0;k<files.size();k++)
	{
	vector <string>* tempVecStr = files[k];
	vector <int> loopNumbers = getLastInts(tempVecStr); //not unsigned so that -1 can be used
	for (unsigned int l=0;l<(*tempVecStr).size();l++)
		{
		if (loopNumbers[l]<firstLoop)
			{
			(*tempVecStr).erase((*tempVecStr).begin()+l);
			}
		}
	}

//reducing to minimum lists with all the same timeNumbers
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

//sorting files		
sort(piFiles.begin(), piFiles.end());
sort(inputsFiles.begin(), inputsFiles.end());

//getting timeNumbers
vector <unsigned long long int> fileNumbers = getInts(piFiles);

//printing filenames
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
	
	//defining the timeNumber
	string timeNumber = currentDateTime();

	//copying a version of mainInputs with timeNumber
	string runInputs = "./data/" + timeNumber + "mainInputs";
	copyFile("mainInputs",runInputs);

	ifstream fin;
	fin.open(inputsFiles[fileLoop].c_str());
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
		V_params = &V1_params;
		dV_params = &dV1_params;
		ddV_params = &ddV1_params;
		epsilon = dE;
		epsilon0 = 0.0;
		}
	else if (pot[0]=='2')
		{
		V_params = &V2_params;
		dV_params = &dV2_params;
		ddV_params = &ddV2_params;
		epsilon = 0.75;
		epsilon0 = 0.7450777428719992;
		}
	else
		{
		cout << "pot option not available, pot = " << pot << endl;
		}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//finding epsilon and root

	//gsl function for dV(phi)
	struct f_gsl_params fparams = { epsilon, A};
	gsl_function_fdf FDF;
	FDF.f = f_gsl;
	FDF.df = df_gsl;
	FDF.fdf = fdf_gsl;
	FDF.params = &fparams;	
	
	//finding roots of dV(phi)=0
	root = minimaFn(&FDF, -3.0, 3.0, 20);
	sort(root.begin(),root.end());
	
	//gsl function for V(root2)-V(root1)-dE
	struct ec_gsl_params ec_params = { A, root[0], root[2], dE};
	gsl_function EC;
	EC.function = &ec_gsl;
	EC.params = &ec_params;
	
	//evaluating epsilon, new root and dE may change slightly
	epsilonFn(&FDF,&EC,&dE,&epsilon,&root);
	
	//evaluating V and a couple of properties of V
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
	comp ergZero = N*a*V(root[0]);
	mass2 = real(ddV(root[0]));
	
	//finding root0 of dV0(phi)=0;
	struct void_gsl_params vparams = {};
	vector<double> root0(3);
	if (pot[0]=='1')
		{
		root0[0] = -1.0; root0[1] = 0.0; root0[2] = 1.0;
		}
	else if (pot[0]=='2')
		{
		gsl_function_fdf DV0DDV0;
		DV0DDV0.f = dV0_gsl;
		DV0DDV0.df = ddV0_gsl;
		DV0DDV0.fdf = dV0ddV0_gsl;
		DV0DDV0.params = &vparams;	
		root0 = minimaFn(&DV0DDV0, -2.0, 2.0, 20);
		sort(root0.begin(),root0.end());
		}
	
	//finding S1
	double S1, S1error;
	gsl_function S1_integrand;
	S1_integrand.function = &s1_gsl;
	S1_integrand.params = &vparams;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
	gsl_integration_qag(&S1_integrand, root0[0], root0[2], DBL_MIN, 1.0e-8, (size_t)1e4, 1e4, w, &S1, &S1error);
	gsl_integration_workspace_free(w);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//other derived quantities
	NT = Na + Nb + Nc;
	R = S1/dE;
	alpha *= R;
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
	
	#pragma omp parallel for	
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
			//system("./negEig"); //negEig now needs timeNumber
			//char * fileNumber = (char *)(fileNumbers[fileLoop]);
			//system(fileNumber); //won't work if inF=='m', as needs fileNumber of pi run
			//cout << "negEig run" << endl;
			cout << "negEigDone==0, run negEig and then set negEigDone=1" << endl;
			}
		negVec = loadVector("data/eigVec.dat",Nb,1);
		ifstream eigFile;
		eigFile.open("data/eigValue.dat");
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
			if ((absolute(theta-minTheta)>DBL_MIN || absolute(Tb-minTb)>DBL_MIN) && loop==0)
				{
				cout << "input Tb: " << Tb << endl;
				cout << "program Tb: " << minTb << endl;
				cout << "input theta: " << theta << endl;
				cout << "program theta: " << minTheta << endl;
				cout << endl;
				}
			if (absolute(maxTheta-minTheta)>DBL_MIN)
				{
				theta = minTheta + (maxTheta - minTheta)*loop/(loops-1.0);
				Gamma = exp(-theta);
				}
			else if (absolute(maxTb-minTb)>DBL_MIN)
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
		double twaction = -pi*epsilon*pow(R,2)/2.0 + pi*R*S1;
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

		//initializing phi (=p)
		vec p(2*N*NT+2);
		if (loop==0)
			{
			p = loadVector(piFiles[fileLoop],NT,2);
			printf("%12s%30s\n","input: ",(piFiles[fileLoop]).c_str());
			}
		else
			{
			unsigned int toLoad = loop-1;
			string loadfile = "./data/" + timeNumber + "mainpi_" + numberToString<unsigned int>(toLoad)+".dat";
			p = loadVector(loadfile,NT,2);
			printf("%12s%30s\n","input: ",loadfile.c_str());
			}
			
		//printing loop name and parameters
		printf("%12s%12s\n","timeNumber: ",timeNumber.c_str());
				
		printParameters();
			
		//very early vector print
		string earlyPrintFile = "data/" + timeNumber + "mainpiE" + numberToString<unsigned int>(loop) + "_" + "0.dat";
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
				char* tempChar = &zmx[zmx.size()-1];
				stringstream ssX(tempChar);
				unsigned int slicesX, slicesT;
				ssX >> slicesX;
				tempChar = &zmt[zmt.size()-1];
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
								chiX(posX+k) = negVec(2*(posCe+k));
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
			if (absolute(normX)<DBL_MIN || absolute(normT)<DBL_MIN)
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
			DDS_to_reserve(0) = N+1;
			DDS_to_reserve(1) = N+1;
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//assigning values to minusDS and DDS and evaluating action
			for (unsigned long int j = 0; j < N*NT; j++)
				{		
				unsigned int t = intCoord(j,0,NT); //coordinates
				unsigned int x = intCoord(j,1,NT);
				comp Dt = DtFn(t);
				comp dt = dtFn(t);
			
				if (absolute(chiX(j))>DBL_MIN) //spatial zero mode lagrange constraint
					{
					DDS.insert(2*j,2*N*NT) = a*chiX(j); 
					DDS.insert(2*N*NT,2*j) = a*chiX(j);
					minusDS(2*j) += -a*chiX(j)*p(2*N*NT);
					minusDS(2*N*NT) += -a*chiX(j)*p(2*j);
					}
					
				if (absolute(chiT(j))>DBL_MIN)
					{
					DDS.coeffRef(2*(j+1),2*N*NT+1) += a*chiT(j); //chiT should be 0 at t=(NT-1) or this line will go wrong
					DDS.coeffRef(2*N*NT+1,2*(j+1)) += a*chiT(j);
					DDS.coeffRef(2*j,2*N*NT+1) += -a*chiT(j);
					DDS.coeffRef(2*N*NT+1,2*j) += -a*chiT(j);
		            minusDS(2*(j+1)) += - a*chiT(j)*p(2*N*NT+1);
		            minusDS(2*j) += a*chiT(j)*p(2*N*NT+1);
		            minusDS(2*N*NT+1) += - a*chiT(j)*(p(2*(j+1))-p(2*j));
					}
					
				if (absolute(theta)<DBL_MIN)
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
					
					if (absolute(theta)<DBL_MIN)
						{
						//DDS.insert(2*j,2*(j+1)+1) = 1.0; //zero imaginary part of time derivative
						//DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
						/////////////////////////////////////equation 1
						for (unsigned int k=0;k<N;k++) //not the simplest b.c.s possible, but continuously related to the theta!=0 ones
							{
							if (absolute(omega(x,k))>DBL_MIN)
								{
								unsigned int m=k*NT;
								DDS.coeffRef(2*j+1,2*m+1) += -2.0*omega(x,k);
								minusDS(2*j+1) += 2.0*omega(x,k)*p(2*m+1);
								}
							}
						//////////////////////////////////////equation 2
						for (unsigned int k=1; k<2*2; k++)
				        	{
				            int sign = pow(-1,k+1);
				            int direc = (int)(k/2.0);
				            if (direc == 0)
				            	{
				                minusDS(2*j) += imag(a/dt)*p(2*j) + real(a/dt)*p(2*j+1);
				                DDS.coeffRef(2*j,2*(j+sign)) += -imag(a/dt);
				                DDS.coeffRef(2*j,2*(j+sign)+1) += -real(a/dt);
				                }
				            else
				            	{
				                unsigned int neighb = neigh(j,direc,sign,NT);
				                minusDS(2*j) += - imag(Dt/a)*p(2*neighb) - real(Dt/a)*p(2*neighb+1);
				                DDS.coeffRef(2*j,2*neighb) += imag(Dt/a);
				                DDS.coeffRef(2*j,2*neighb+1) += real(Dt/a);
				                }
				            }
				        comp temp0 = a/dt - 2.0*Dt/a;
				        double temp1 = imag(a*Dt*dV(Cp(j)));
				        double temp2 = - 3.0*real(Dt)*a*p(2*j)*p(2*j+1);
				        
				        minusDS(2*j) += -imag(temp0)*p(2*j) - real(temp0)*p(2*j+1) + imag(a*Dt*dV(Cp(j)));
				        //////////////////////UNFINISHED/////////////////////////////
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
		            	comp temp1 = a*Dt*(2.0*Cp(j)/pow(a,2.0) + dV(Cp(j)) + dVr(Cp(j)));//dV terms should be small as linearised

				        minusDS(2*j) += real(temp1 - temp0*Cp(j));
				        minusDS(2*j+1) += imag(temp1 - temp0*Cp(j));

                        for (unsigned int k=0; k<N; k++)
                        	{
                        	//dropped two terms here, as had imag(real()) or vice versa
							unsigned int m = x*NT; //equals j
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
                    comp dtd = dt;
                    if (sign==-1) {dtd = dtFn(t-1);}
                    if (direc == 0)
                    	{
                        minusDS(2*j) += real(a*Cp(j+sign)/dtd);
                        minusDS(2*j+1) += imag(a*Cp(j+sign)/dtd);
                        DDS.insert(2*j,2*(j+sign)) = -real(a/dtd);
                        DDS.insert(2*j,2*(j+sign)+1) = imag(a/dtd);
                        DDS.insert(2*j+1,2*(j+sign)) = -imag(a/dtd);
                        DDS.insert(2*j+1,2*(j+sign)+1) = -real(a/dtd);
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
		            comp temp0 = a*(1.0/dt + 1.0/dtFn(t-1));
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
			double momTest = E*b/Num*pi; //perhaps should have a not b here
			mom_test.push_back(momTest);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			string prefix = "./data/" + timeNumber;
			string suffix = numberToString<unsigned int>(loop)+"_"+numberToString<unsigned int>(runs_count)+".dat";
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
				string piEarlyFile = prefix + "mainpiE"+suffix;
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
				cout << endl;
				cout << "ergTest = " << erg_test.back() << endl;
				}
		
		//checking lattice small enough
		if (mom_test.back()>closenessP)
			{
			cout << "momTest = "<< mom_test.back()  << endl;
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
		
		string prefix = "./data/" + timeNumber;
		string suffix = "_" + numberToString<unsigned int>(fileLoop)+"_" + numberToString<unsigned int>(fileLoop)+ ".dat";
	
		//copying a version of inputs with timeNumber and theta changed
		string runInputs = prefix + "inputsM"+ numberToString<unsigned int>(fileLoop) + "_" + numberToString<unsigned int>(loop); //different suffix
		if (absolute(maxTheta-minTheta)>DBL_MIN)
			{
			changeInputs(runInputs, "theta", numberToString<double>(theta));
			}
		else if (absolute(maxTb-minTb)>DBL_MIN)
			{
			changeInputs(runInputs, "Tb", numberToString<double>(Tb));
			}
		else
			{
			copyFile("inputs",runInputs);
			}
	
		//printing output phi
		string tpifile =  prefix + "mainpi"+suffix;
		printVector(tpifile,p);
		gp(tpifile,"repi.gp");
	
		//printing output minusDS				
		//string minusDSfile = prefix + "mainminusDS"+suffix;
		//printVector(minusDSfile,minusDS);
				
		//printing output DDS
		//string DDSfile = prefix + "mainDDS"+suffix;
		//printSpmat(DDSfile,DDS);
	
		//printing linNumVec
		//string linNumFile = prefix + "mainlinNum"+suffix;
		//linNum.conservativeResize(Na);
		//simplePrintVector(linNumFile,linNum);
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
		
		} //ending parameter loop
	} //ending file loop

return 0;
}

