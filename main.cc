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

//getting list of relevant data files
system("dir ./data/* > dataFiles");
vector<string> filenames, piFiles, inputsFiles, eigenvectorFiles, eigenvalueFiles;
filenames = readDataFiles(minFile,maxFile);
piFiles = findStrings(filenames,"tpip");
inputsFiles = findStrings(filenames,"inputs");
eigenvectorFiles = findStrings(filenames,"eigVec");
eigenvalueFiles = findStrings(filenames,"eigValue");
if (piFiles.size()!=inputsFiles.size() || piFiles.size()!=eigenvectorFiles.size() || piFiles.size()!=eigenvalueFiles.size())
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
		negVec = loadVector(eigenvectorFiles[fileLoop],Nb);
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
		
	//defining the timeNumber
	string timeNumber = currentDateTime();

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
	
		//defining the energy vector
		cVec erg(NT-1);
		
		//defining the action
		comp action = twaction;
		
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

		//initializing phi (=p)
		vec p(2*N*NT+1);
		if (loop==0)
			{
			p = loadVector(piFiles[fileLoop],NT);
			}
		else
			{
			unsigned int toLoad = loop-1;
			string loadfile = "./data/" + timeNumber + "main" + to_string(toLoad)+".dat";
			p = loadVector(loadfile,NT);
			}
			
		//very early vector print
		string earlyPrintFile = "data/" + timeNumber + "mainE"+ to_string(loop) + "0.dat";
		printVectorB(earlyPrintFile,p);
		
		//defining complexified vector Cp
		cVec Cp(NT*N);
		Cp = vecComplex(p,N*NT);
	
		//defining DDS and minusDS
		spMat DDS(2*N*NT+1,2*N*NT+1);
		vec minusDS(2*N*NT+1);
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//beginning newton-raphson loop	
		while ((sol_test.back()>closenessS || solM_test.back()>closenessSM || runs_count<min_runs))
			{
			runs_count++;
			
			}
		}
	}

return 0;
}

