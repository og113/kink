#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include "pf.h"
#include "files.h"

#define MIN_NUMBER 1.0e-16

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
if (fmainin.is_open()) {
	string line;
	unsigned int lineNumber = 0;
	while(getline(fmainin,line)) {
		if(line[0] == '#') continue;
		else if(line.empty()) continue;
		istringstream ss(line);
		if (lineNumber==0) {
			ss >> inF >> minFile >> maxFile >> firstLoop >> lastLoop;
			lineNumber++;
		}
		else if (lineNumber==1) {
			ss >> minTb >> maxTb >> minTheta >> maxTheta >> loops >> zmt >> zmx >> bds;
			lineNumber++;
		}
	}
}
else {
	cerr << "unable to open mainInputs" << endl;
	return 1;
}
fmainin.close();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//sorting out files to load

//getting list of relevant data files, with timeNumbers between minFile and maxFile
int systemCall = system("dir ./data/* > dataFiles");
if (systemCall==-1){
	cerr << "system call failure, finding dataFiles" << endl;
	return 1;
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
		for (unsigned int j=0; j<eigenvalueFiles.size(); j++) cout << eigenvalueFiles[j] << endl;
		cout << "eigenvectorFiles.size() = " << eigenvectorFiles.size() << endl;
		for (unsigned int j=0; j<eigenvectorFiles.size(); j++) cout << eigenvectorFiles[j] << endl;
		}
	}
else if (inF.compare("m")==0)
	{
	piFiles = findStrings(filenames,"mainpi_");
	inputsFiles = findStrings(filenames,"inputsM");
	}
else
	{
	cerr << "inF error" << endl;
	return 1;
	}

//picking subset based on loopNumbers being above firstLoop
vector <vector<string>*> files;
files.reserve(2);
files.push_back(&piFiles); files.push_back(&inputsFiles);
if (files.size()==0)
	{
	cerr << "no files found" << endl;
	return 1;
	}
for (unsigned int k=0;k<files.size();k++)
	{
	vector <string>* tempVecStr = files[k];
	vector <int> loopNumbers = getLastInts(tempVecStr); //not unsigned so that -1 can be used
	unsigned int l=0;
	while(l<(*tempVecStr).size())
		{
		if (loopNumbers[l]<firstLoop || loopNumbers[l]>lastLoop)
			{
			(*tempVecStr).erase((*tempVecStr).begin()+l);
			loopNumbers.erase(loopNumbers.begin()+l);
			}
		else l++;
		}
	}

//reducing to minimum lists with all the same timeNumbers
if (piFiles.size()!=inputsFiles.size() || piFiles.size()==0)
	{
	for (unsigned int j=0;j<files.size();j++)
		{
		for (unsigned int k=0;k<files.size();k++) *files[j] = reduceTo(*files[j],*files[k]);
		}
	if (piFiles.size()!=inputsFiles.size() || piFiles.size()==0)
		{
		cout << "required files not available" << endl;
		cout << "piFiles.size() = " << piFiles.size() << endl;
		cout << "inputsFiles.size() = " << inputsFiles.size() << endl;
		return 1;
		}
	}

//sorting files		
sort(piFiles.begin(), piFiles.end());
sort(inputsFiles.begin(), inputsFiles.end());

//getting timeNumbers
vector <unsigned long long int> fileNumbers = getInts(piFiles);

//printing filenames
cout << endl;
for (unsigned int j=0; j<inputsFiles.size();j++) cout << inputsFiles[j] << " " << piFiles[j] << endl;
cout << endl;

//defining the timeNumber
string timeNumber = currentDateTime();
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//beginning file loop
for (unsigned int fileLoop=0; fileLoop<piFiles.size(); fileLoop++)
	{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//getting variables from inputs

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
			if(line[0] == '#') continue;
			else if(line.empty()) continue;
			istringstream ss(line);
			if (lineNumber==0)
				{
				ss >> N >> Na >> Nb >> Nc >> dE >> LoR >> Tb >> theta;
				lineNumber++;
				}
			else if (lineNumber==1)
				{
				ss >> aq.inputChoice >> aq.inputTimeNumber >> aq.inputLoop >> aq.totalLoops >> aq.loopChoice >> aq.minValue >> aq.maxValue >> aq.printChoice >> aq.printRun;
				lineNumber++;
				}
			else if(lineNumber==2)
				{
				ss >> alpha >> open >> amp >> pot >> A >> reg >> cutoff;
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
		cerr << "unable to open " << inputsFiles[fileLoop] << endl;
		return 1;
		}
	fin.close();
	inP = aq.inputChoice;

	string loop_choice = aq.loopChoice; //just so that we don't have two full stops when comparing strings
	string print_choice = aq.printChoice;
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//potential function
	
//potential functions
	if (pot[0]=='1')
		{
		neigh = &periodic;
		simpleSpace = &simpleSpaceBox;
		V = &V1c;
		dV = &dV1c;
		ddV = &ddV1c;
		Vd = &V1;
		dVd = &dV1;
		ddVd = &ddV1;
		epsilon0 = 0.0;
		epsilon = dE;
		r0 = 0.0;
		}
	else if (pot[0]=='2')
		{
		neigh = &periodic;
		simpleSpace = &simpleSpaceBox;
		V = &V2c;
		dV = &dV2c;
		ddV = &ddV2c;
		Vd = &V2;
		dVd = &dV2;
		ddVd = &ddV2;
		epsilon0 = 0.74507774287199924;
		epsilon = 0.75;
		r0 = 0.0;
		}
	else if (pot[0]=='3')
		{
		neigh = &spherical;
		simpleSpace = &simpleSpaceSphere;
		V = &V3c;
		dV = &dV3c;
		ddV = &ddV3c;
		Vd = &V3;
		dVd = &dV3;
		ddVd = &ddV3;
		epsilon0 = 0.0;
		epsilon = 0.0;
		r0 = MIN_NUMBER;
		}
	else
		{
		cout << "pot option not available, pot = " << pot << endl;
		}
	paramsV  = {epsilon, A};
	paramsV0 = {epsilon0, A};
	paramsVoid = {};
	
if (cutoff>N || cutoff<1) {
	cerr << "cutoff(" << cutoff << ") too small" << endl;
	return 1;
}

		
//finding epsilon and root
vector<double> minima0(2);
	
	if (pot[0]!='3')
		{
		//gsl function for dV(phi)
		gsl_function F;
		F.function = Vd;
		F.params = &paramsV;	
	
		//finding preliminary roots of dV(phi)=0
		minima[0] = brentMinimum(&F, -1.0, -3.0, 0.0);
		minima[1] = brentMinimum(&F, 1.2, 0.5, 3.0);
	
		//gsl function for V(root2)-V(root1)-dE
		struct ec_params ec_params = { A, minima[0], minima[1], dE};
		gsl_function EC;
		EC.function = &ec;
		EC.params = &ec_params;

		//evaluating epsilon, new root and dE may change slightly
		epsilonFn(&F,&EC,&dE,&epsilon,&minima);
	
		//evaluating mass about false vac
		mass2 = ddVd(minima[0],&paramsV);
	
		//finding root0 of dV0(phi)=0;
		vector<double> minima0(3);
		if (pot[0]=='1')
			{
			minima0[0] = -1.0; minima0[1] = 1.0;
			}
		else if (pot[0]=='2')
			{
			gsl_function V0;
			V0.function = Vd;
			V0.params = &paramsV0;	
			minima0[0] = brentMinimum(&V0, -1.0, -3.0, 0.0);
			minima0[0] = brentMinimum(&V0, 1.2, 0.5, 3.0);
			struct ec_params ec0_params = { A, minima0[0], minima0[1], 0.0};
			gsl_function EC0;
			EC0.function = &ec;
			EC0.params = &ec0_params;
			double dE0 = 0.0;
			epsilonFn(&V0,&EC0,&dE0,&epsilon0,&minima0);
			}
	
		//finding S1
		double S1error;
		gsl_function S1_integrand;
		S1_integrand.function = &s1Integrand;
		S1_integrand.params = &paramsV0;
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
		gsl_integration_qag(&S1_integrand, minima0[0], minima0[1], MIN_NUMBER, 1.0e-8, 1e4, 4, w, &S1, &S1error);
		gsl_integration_workspace_free(w);
		if (S1error>1.0e-8) { cout << "S1 error = " << S1error << endl;}
		}
	else
		{
		mass2 = 1.0;
		minima[0] = 0.0; //only one minimum
		R = 10.0; alpha = 0.0; //not used
		}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//other derived quantities
	NT = Na + Nb + Nc;
	double twaction;
	L = LoR*R;
	if (pot[0]!='3')
		{
		R = S1/dE;
		twaction = -pi*epsilon*pow(R,2)/2.0 + pi*R*S1;
		alpha *= R;
		L = LoR*R;
		if (Tb<R)
			{
			angle = asin(Tb/R);
			double Ltemp = 1.5*(1.5*Tb*tan(angle));
			if (Ltemp<L) L=Ltemp;//making sure to use the smaller of the two possible Ls
			}
		}
	else
		{
		twaction = 8.0*pow(pi,2.0)/3.0; //massless limit
		//no need for R or alpha
		}
	Gamma = exp(-theta);
	vec negVec(2*N*Nb+1);
	a = (L-r0)/(N-1.0);
	b = Tb/(Nb-1.0);
	if (a>0.5*pow(mass2,0.5) || b>0.5*pow(mass2,0.5)) {cerr << endl << "a = " << a << " , b = " << b << endl << endl;}
	Ta = b*Na;
	Tc = b*Nc;
	double ergZero = 0.0;
	if (pot[0]!='3') ergZero = N*a*Vd(minima[0],&paramsV);
	
	//determining number of runs
	closenessA = 1.0;
	closenessS = 1.0e-6;
	closenessSM = 1.0e-5;
	closenessD = 0.7;
	closenessC = 1.0e-16*N*NT;
	closenessCon = 1.0/3.0;
	closenessL = 1.0e-2;
	closenessT = 5.0e-2;
	closenessP = 0.5;
	closenessR = 1.0e-2;
	closenessIE = 1.0e-5;
	closenessCL = 5.0e-2;
	closenessON = 5.0e-2;
	closenessAB = 1.0e-2;
	closenessLR = 1.0e-12;
	closenessABNE = 5.0e-2;
	
	//lambda functions for pot_r
auto Vr = [&] (const comp & phi)
	{
	return -ii*reg*VrFn(phi,minima[0],minima[1]);
	};
auto dVr = [&] (const comp & phi)
	{
	return -ii*reg*dVrFn(phi,minima[0],minima[1]);
	};
auto ddVr = [&] (const comp & phi)
	{
	return -ii*reg*ddVrFn(phi,minima[0],minima[1]);
	};

	//deterimining omega matrices for fourier transforms in spatial direction
	mat omega_m1(N,N), 	omega_m1_cutoff(N,N);
	mat omega_0(N,N), 	omega_0_cutoff(N,N);
	mat omega_1(N,N), 	omega_1_cutoff(N,N);
	mat omega_2(N,N), 	omega_2_cutoff(N,N);
	omega_m1 = Eigen::MatrixXd::Zero(N,N);
	omega_0 = Eigen::MatrixXd::Zero(N,N);
	omega_1 = Eigen::MatrixXd::Zero(N,N);
	omega_2 = Eigen::MatrixXd::Zero(N,N);
	omega_m1_cutoff = Eigen::MatrixXd::Zero(N,N);
	omega_0_cutoff = Eigen::MatrixXd::Zero(N,N);
	omega_1_cutoff = Eigen::MatrixXd::Zero(N,N);
	omega_2_cutoff = Eigen::MatrixXd::Zero(N,N);
	vec eigenValues(N);
	mat eigenVectors(N,N); //eigenvectors correspond to columns of this matrix - so the nth eigenvector is v_n(j) = eigenVectors(j,n)
	bool approxOmega = false;
	if (!approxOmega) {
		mat h(N,N);
		h = hFn(N,a,mass2);
		Eigen::SelfAdjointEigenSolver<mat> eigensolver(h);
		if (eigensolver.info() != Eigen::Success) {
			cerr << "h eigensolver failed" << endl;
			cerr << "N = " << N << ", a = " << a << ", mass2 = " << mass2 << endl;
		}
		else {
			eigenValues = eigensolver.eigenvalues();
			eigenVectors = eigensolver.eigenvectors(); //automatically normalised to have unit norm
		}
	}
	else {
		double normalisation = sqrt(2.0/(N-1.0));
		for (unsigned int l=0; l<N; l++) {
			eigenValues(l) = 1.0+pow(2.0*sin(pi*l/(N-1.0)/2.0)/a,2.0);
			for (unsigned int m=0; m<N; m++) {
				if (pot[0]=='3') eigenVectors(l,m) = normalisation*sin(pi*l*m/(N-1.0));
				else			 eigenVectors(l,m) = normalisation*cos(pi*l*m/(N-1.0));
			}
		}
	}
	double djdk;	
	for (unsigned int j=0; j<N; j++) {
		for (unsigned int k=0; k<N; k++) {
			for (unsigned int l=0; l<N; l++) {
				if (approxOmega) djdk = a;
				else 			 djdk = sqrt(DxFn(j)*DxFn(k));
				omega_m1(j,k) += djdk*pow(eigenValues(l),-0.5)*eigenVectors(j,l)*eigenVectors(k,l);
				omega_0(j,k)  += djdk*eigenVectors(j,l)*eigenVectors(k,l);
				omega_1(j,k)  += djdk*pow(eigenValues(l),0.5)*eigenVectors(j,l)*eigenVectors(k,l);
				omega_2(j,k)  += djdk*eigenValues(l)*eigenVectors(j,l)*eigenVectors(k,l);
				if (l<cutoff) {
					omega_m1_cutoff(j,k) += djdk*pow(eigenValues(l),-0.5)*eigenVectors(j,l)*eigenVectors(k,l);
					omega_0_cutoff(j,k)  += djdk*eigenVectors(j,l)*eigenVectors(k,l);
					omega_1_cutoff(j,k)	 += djdk*pow(eigenValues(l),0.5)*eigenVectors(j,l)*eigenVectors(k,l);
					omega_2_cutoff(j,k)	 += djdk*eigenValues(l)*eigenVectors(j,l)*eigenVectors(k,l);
				}
			}
		}
	}

	if (zmt[0]=='n' || zmx[0]=='n')
		{
		if (pot[0]=='3') {
			vec negVecFull = loadVectorColumn("data/stable/sphaleronEigVec.dat",0); // nb may need to check that r1 is correct
			negVec = interpolate1d(negVecFull,negVecFull.size(),N);
		}
		else
			{
			if (negEigDone==0)
				{
				//system("./negEig"); //negEig now needs timeNumber
				//char * fileNumber = (char *)(fileNumbers[fileLoop]);
				//system(fileNumber); //won't work if inF=='m', as needs fileNumber of pi run
				//cout << "negEig run" << endl;
				cerr << "negEigDone==0, run negEig and then set negEigDone=1" << endl;
				return 1;
				}
		string tempPot = pot.substr(0,1);
		string eigVecFilename = "./data/stable/eigVec_" + tempPot + ".dat";
		unsigned int fileLength = countLines(eigVecFilename);
		if (fileLength==(N*Nb+1)) negVec = loadVector(eigVecFilename,Nb,N,1);
		else if (fileLength % 2) //if its odd
			{
			string eigVecInputs = "data/stable/eigVecInputs" + tempPot + ".dat";
			unsigned int NEig, NtEig;
			ifstream fin;
			fin.open(eigVecInputs.c_str());
			if (!fin.good()) {
				cerr << eigVecInputs << " didn't open successfully" << endl;
				return 1;
			}
			if (fin.is_open())
				{
				string line, tempStr;
				while(getline(fin,line))
					{
					if(line[0] == '#') continue;
					if(line.empty()) continue;
					istringstream ss(line);
					ss >> NEig >> tempStr >> NtEig;
					break;
					}
				}
			else {
				cerr << "unable to open " << eigVecInputs << endl;
				return 1;
			}
			fin.close();
			vec tempEig = loadVector(eigVecFilename,NtEig,NEig,1);
			negVec = interpolate(tempEig,NtEig,NEig,Nb,N);
			}
		else
			{
			cerr << "eigVec not the right length, cannot interpolate" << endl;
			return 1;
			}
			ifstream eigFile;
			eigFile.open("data/stable/eigValue.dat");
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//beginning theta or Tb loop
	for (unsigned int loop=0; loop<loops; loop++)
		{
		if (loops>1)
			{
			if ((abs(theta-minTheta)>MIN_NUMBER || abs(Tb-minTb)>MIN_NUMBER) && loop==0)
				{
				cout << "input Tb     : " << Tb << endl;
				cout << "program Tb   : " << minTb << endl;
				cout << "input theta  : " << theta << endl;
				cout << "program theta: " << minTheta << endl;
				cout << endl;
				}
			if (abs(maxTheta-minTheta)>MIN_NUMBER)
				{
				theta = minTheta + (maxTheta - minTheta)*loop/(loops-1.0);
				Gamma = exp(-theta);
				}
			else if (abs(maxTb-minTb)>MIN_NUMBER)
				{
				Tb = minTb + (maxTb - minTb)*loop/(loops-1.0);
				changeDouble ("Tb",Tb);
				}
			else
				{
				cout << "nothing to loop over, set loops to 1" << endl;
				loops = 1;
				}
			}
		
		//defining a time and starting the clock
		clock_t time;
		time = clock();
	
		//defining energy and number vectors
		cVec erg(NT);
		cVec linErg(NT);
		cVec linNum(NT);
		cVec linErgOffShell(NT);
		cVec linNumOffShell(NT);
		cVec derivErg(NT), potErg(NT);
		comp linErgContm, linNumContm;
		comp linNumAB, linErgAB;
		
		//defining the action and bound and W
		comp action = ii*twaction;
		double bound;
		double W;
		double E;
		double E_exact;
		double Num;
		
		//defining some quantities used to stop the Newton-Raphson loop when action stops varying
		comp action_last = action;
		unsigned int runs_count = 0;
		unsigned int min_runs = 3;
		vector<double> action_test(1);	action_test[0] = 1.0;
		vector<double> sol_test(1);	sol_test[0] = 1.0;
		vector<double> solM_test(1);	solM_test[0] = 1.0;
		vector<double> delta_test(1); delta_test[0] = 1.0;
		vector<double> calc_test(1); calc_test[0] = 1.0;
	 	vector<double> conserv_test(1); conserv_test[0] = 1.0;
		vector<double> lin_test(1); lin_test[0] = 1.0;
		vector<double> true_test(1); true_test[0] = 1.0;
		vector<double> mom_test(1); mom_test[0] = 1.0;
		vector<double> reg_test(1); reg_test[0] = 1.0;
		vector<double> imErg_test(1); imErg_test[0] = 1.0;
		vector<double> onShell_test(1); onShell_test[0] = 1.0;
		vector<double> AB_test(1); AB_test[0] = 1.0;
		vector<double> linRep_test(1); linRep_test[0] = 1.0;
		vector<double> ABNE_test(1); ABNE_test[0] = 1.0;

		//initializing phi (=p)
		vec p(2*N*NT+2);
		if (loop==0)
			{
			p = loadVector(piFiles[fileLoop],NT,N,2);
			printf("%12s%30s\n","input: ",(piFiles[fileLoop]).c_str());
			}
		else
			{
			string loadfile = "./data/" + timeNumber + "mainpi_" + numberToString<unsigned int>(fileLoop) + "_" + numberToString<unsigned int>(loop-1)+".dat";
			p = loadVector(loadfile,NT,N,2);
			printf("%12s%30s\n","input: ",loadfile.c_str());
			}
			
		//printing loop name and parameters
		printf("%12s%12s\n","timeNumber: ",timeNumber.c_str());
				
		printParameters();
		//printMoreParameters();
			
		//very early vector print
		string earlyPrintFile = "data/" + timeNumber + "mainpiE" + "_" + numberToString<unsigned int>(fileLoop) + "_" + numberToString<unsigned int>(loop) + "_" + "0.dat";
		printVector(earlyPrintFile,p);
		
		//defining complexified vector Cp
		cVec Cp(NT*N);
		Cp = vecComplex(p,NT*N);
	
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
				unsigned int slicesX, slicesT;
				if (getLastInt(zmx)<0) {
					cerr << "getLastInt error with zmx = " << zmx << endl;
					return 1;
				}
				if (getLastInt(zmt)<0) {
					cerr << "getLastInt error with zmt = " << zmt << endl;
					return 1;
				}
				slicesX = getLastInt(zmx);
				slicesT = getLastInt(zmt);
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
							if (zmt[0]=='n' && pot[0]!='3') 	chiT(posT+k) = negVec(2*(posCe+k));
							else if (zmt[0]=='n' && pot[0]=='3')
								{
								double r = r0 + j*a;
								chiT(posT+k) = negVec(j)*r;
								}
							else if (zmt[0]=='d')				chiT(posT+k) = p(2*(posT+k+1))-p(2*(posT+k));
							else
								{
								cerr << "choice of zmt not allowed" << endl;
								return 1;
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
							if (zmx[0]=='n' && pot[0]!='3')		chiX(posX+k) = negVec(2*(posCe+k));
							else if (zmx[0]=='n' && pot[0]=='3')
								{
								double r = r0 + j*a;
								chiX(posX+k) = negVec(j)*r;
								}
							else if (zmx[0]=='d' && pot[0]!='3') chiX(posX+k) = p(2*neigh(posX+k,1,1,NT,N))-p(2*neigh(posX+k,1,-1,NT,N));
							else
								{
								cerr << "choice of zmx not allowed" << endl;
								return 1;
								}
							}
						}
					}
				}
			double normX = chiX.norm();
			double normT = chiT.norm();
			normT = pow(normT,0.5);
			if (abs(normX)<MIN_NUMBER || abs(normT)<MIN_NUMBER)
				{
				cerr << "norm of chiX = " << normX << ", norm of chiT = " << normT << endl;
				}
			chiX = chiX/normX;
			chiT = chiT/normT;
			
			// allocating memory for DS, DDS
			minusDS = Eigen::VectorXd::Zero(2*N*NT+2); //initializing to zero
			DDS.setZero(); //just making sure
			Eigen::VectorXi DDS_to_reserve(2*N*NT+2);//number of non-zero elements per column
			DDS_to_reserve = Eigen::VectorXi::Constant(2*N*NT+2,13);
			if (abs(theta)<2.0e-16)
				{
				DDS_to_reserve(0) = N+3;
				DDS_to_reserve(1) = 11;
				}
			else
				{
				DDS_to_reserve(0) = N+11;
				DDS_to_reserve(1) = N+11;
				}
			DDS_to_reserve(2*N*NT-2) = 4;
			DDS_to_reserve(2*N*NT-1) = 4;
			DDS_to_reserve(2*N*NT) = N*(zmx.size()-2);
			DDS_to_reserve(2*N*NT+1) = 2*N*(zmt.size()-2);
			DDS.reserve(DDS_to_reserve);
			
			//initializing to zero
			comp kineticS = 0.0;
			comp kineticT = 0.0;
			comp potV = 0.0;
			comp pot_r = 0.0;
			bound = 0.0;
			erg = Eigen::VectorXcd::Constant(NT,-ergZero);
			linErg = Eigen::VectorXcd::Zero(NT);
			linErgOffShell = Eigen::VectorXcd::Zero(NT);
			linNum = Eigen::VectorXcd::Zero(NT);
			linNumOffShell = Eigen::VectorXcd::Zero(NT);
			derivErg = Eigen::VectorXcd::Zero(NT);
			potErg = Eigen::VectorXcd::Constant(NT,-ergZero);
			linErgContm = 0.0;
			linNumContm = 0.0;
			
			//testing that the potential term is working for pot3
			if (pot[0]=='3' && false)
				{
				comp Vtrial = 0.0, Vcontrol = 0.0;
				for (unsigned int j=0; j<N; j++)
					{
					double r = r0 + j*a;
					paramsV  = {r, 0.0};
					Vcontrol += pow(p(2*j*Nb),2.0)/2.0 - pow(p(2*j*Nb),4.0)/4.0/pow(r,2.0);
					Vtrial += V(p(2*j*Nb));
					}
				double potTest = pow(pow(real(Vcontrol-Vtrial),2.0) + pow(imag(Vcontrol-Vtrial),2.0),0.5);
				cout << "potTest = " << potTest << endl;
				}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//assigning values to minusDS and DDS and evaluating action
			//#pragma omp parallel for
			for (unsigned long int j = 0; j < N*NT; j++)
				{		
				unsigned int t 			= intCoord(j,0,NT); //coordinates
				unsigned int x 			= intCoord(j,1,NT);
				unsigned int neighPosX 	= neigh(j,1,1,NT,N);
				
				comp Dt 		= DtFn(t);
				comp dt 		= dtFn(t);
				comp dtm 		= (t>0? dtFn(t-1): b);
				double Dx 		= DxFn(x);
				double dx 		= dxFn(x);
				double dxm 		= (x>0? dxFn(x-1): a);
				
				if (pot[0]=='3') paramsV  = {r0+x*a, 0.0};
			
				if (abs(chiX(j))>MIN_NUMBER && pot[0]!='3') //spatial zero mode lagrange constraint
					{
					DDS.insert(2*j,2*N*NT) 			= Dx*chiX(j); 
					DDS.insert(2*N*NT,2*j) 			= Dx*chiX(j);
					minusDS(2*j) 					+= -Dx*chiX(j)*p(2*N*NT);
					minusDS(2*N*NT) 				+= -Dx*chiX(j)*p(2*j);
					}
					
				if (abs(chiT(j))>MIN_NUMBER && t<(NT-1))
					{
					DDS.coeffRef(2*(j+1),2*N*NT+1) 	+= Dx*chiT(j); //chiT should be 0 at t=(NT-1) or this line will go wrong
					DDS.coeffRef(2*N*NT+1,2*(j+1)) 	+= Dx*chiT(j);
					DDS.coeffRef(2*j,2*N*NT+1) 		+= -Dx*chiT(j);
					DDS.coeffRef(2*N*NT+1,2*j) 		+= -Dx*chiT(j);
		            minusDS(2*(j+1)) 				+= -Dx*chiT(j)*p(2*N*NT+1);
		            minusDS(2*j) 					+= Dx*chiT(j)*p(2*N*NT+1);
		            minusDS(2*N*NT+1) 				+= -Dx*chiT(j)*(p(2*(j+1))-p(2*j));
					}
					
				if (t<(NT-1))
					{
					for (unsigned int k=0;k<N;k++)
						{
						unsigned int l = k*NT+t;
						linErgOffShell(t) += 0.5*( omega_2_cutoff(x,k)*( p(2*l)-minima[0] )*( p(2*j)-minima[0] ) \
										+ omega_0_cutoff(x,k)*( p(2*(l+1))-p(2*l) )*( p(2*(j+1))-p(2*j) )/pow(dt,2.0));
						linNumOffShell(t) += 0.5*(omega_1_cutoff(x,k)*( p(2*l)-minima[0] )*( p(2*j)-minima[0] ) \
										+ omega_m1_cutoff(x,k)*( p(2*(l+1))-p(2*l) )*( p(2*(j+1))-p(2*j) )/pow(dt,2.0));
						}
					}
					
				if (abs(theta)>MIN_NUMBER)
					{
					for (unsigned int k=0;k<N;k++)
						{
						unsigned int l = k*NT+t;
						linNum(t) += 2.0*Gamma*omega_1_cutoff(x,k)*( (p(2*l)-minima[0])*(p(2*j)-minima[0])/pow(1.0+Gamma,2.0)\
								 + p(2*j+1)*p(2*l+1)/pow(1.0-Gamma,2.0));
						linErg(t) += 2.0*Gamma*omega_2_cutoff(x,k)*( (p(2*l)-minima[0])*(p(2*j)-minima[0])/pow(1.0+Gamma,2.0)\
								+ p(2*j+1)*p(2*l+1)/pow(1.0-Gamma,2.0));
						}
					}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//boundaries
				if (pot[0]=='3' && x==(N-1))
					{
					DDS.insert(2*j,2*j) 	= 1.0; // p=0 at r=R
					DDS.insert(2*j+1,2*j+1) = 1.0;
					}
				else if (pot[0]=='3' && x==0)
					{
					kineticS 				+= 	Dt*pow(Cp(neighPosX),2.0)/dx/2.0;
					derivErg(t) 			+= 	pow(Cp(neighPosX),2.0)/dx/2.0; //n.b. the Dt/dt difference is ignored for erg(t)
					erg(t) 					+= 	pow(Cp(neighPosX),2.0)/dx/2.0;
					
					DDS.insert(2*j,2*j) 	= 1.0; // p=0 at r=0
					DDS.insert(2*j+1,2*j+1) = 1.0;
					}			
				else if (t==(NT-1))
					{
					kineticS 				+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
					potV 					+= Dt*Dx*V(Cp(j));
					pot_r 					+= Dt*Dx*Vr(Cp(j));
					erg(t) 					+= pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0 + Dx*V(Cp(j)) + Dx*Vr(Cp(j));
					derivErg(t) 			+= pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
					potErg(t) 				+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
				
					DDS.insert(2*j,2*(j-1)+1) = 1.0; //zero imaginary part of time derivative
					DDS.insert(2*j+1,2*j+1)   = 1.0; //zero imaginary part
					}
				else if (t==0)
					{
					kineticT 	+= Dx*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					derivErg(t) += Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0;
					erg(t) 		+= Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0;
					
					if (bds.compare("ft")==0) // this actually requires that p(0,x)=p(0,k) (t=0)
						{
						if (abs(theta)<MIN_NUMBER)
							{
							for (unsigned int k=0; k<N; k++)
								{
								unsigned int m = k*NT;
								minusDS(2*j)				+= 2.0*sqrt(Dx)*eigenVectors(x,k)*eigenValues(k)*p(2*m+1);
								minusDS(2*j+1)				+= sqrt(Dx)*eigenVectors(x,k)*(p(2*(m+1)+1)-p(2*m+1))/real(dt);
								DDS.coeffRef(2*j,2*m+1)		+= -2.0*sqrt(Dx)*eigenVectors(x,k)*eigenValues(k);
								DDS.coeffRef(2*j+1,2*(m+1)+1)+= -sqrt(Dx)*eigenVectors(x,k)/real(dt);
								DDS.coeffRef(2*j+1,2*m+1) 	+= sqrt(Dx)*eigenVectors(x,k)/real(dt);
								}
							}
						else
							{
							for (unsigned int k=0; k<N; k++)
								{
								unsigned int m = k*NT;
								minusDS(2*j)				+= theta*sqrt(Dx)*eigenVectors(x,k)*((p(2*(m+1))-p(2*m))/real(dt) \
																+ ((1.0+Gamma)/(1.0-Gamma))*eigenValues(k)*p(2*m+1));
								minusDS(2*j+1)				+= sqrt(Dx)*eigenVectors(x,k)*((p(2*(m+1)+1)-p(2*m+1))/real(dt) \
																- ((1.0-Gamma)/(1.0+Gamma))*eigenValues(k)*p(2*m));
								DDS.coeffRef(2*j,2*m)		+= theta*sqrt(Dx)*eigenVectors(x,k)/real(dt);
								DDS.coeffRef(2*j,2*(m+1))	+= -theta*sqrt(Dx)*eigenVectors(x,k)/real(dt);
								DDS.coeffRef(2*j,2*m+1)		+= -theta*sqrt(Dx)*eigenVectors(x,k)*eigenValues(k)*(1.0+Gamma)/(1.0-Gamma);
								DDS.coeffRef(2*j+1,2*m+1)	+= sqrt(Dx)*eigenVectors(x,k)/real(dt);
								DDS.coeffRef(2*j+1,2*(m+1)+1)+= -sqrt(Dx)*eigenVectors(x,k)/real(dt);
								DDS.coeffRef(2*j+1,2*m)		+= +sqrt(Dx)*eigenVectors(x,k)*eigenValues(k)*(1.0-Gamma)/(1.0+Gamma);
								}
							}
						}
					else if (bds.compare("jd")==0)
						{
						minusDS(2*j+1) 					+= Dx*imag(Cp(j+1)/dt);
					    DDS.coeffRef(2*j+1,2*(j+1)) 	+= -imag(Dx/dt);
					    DDS.coeffRef(2*j+1,2*(j+1)+1)	+= -real(Dx/dt);
					    minusDS(2*j+1) 					+= imag(-Dx*Cp(j)/dt);
					    DDS.coeffRef(2*j+1,2*j) 		+= imag(Dx/dt);
					    DDS.coeffRef(2*j+1,2*j+1) 		+= real(Dx/dt);
					    if (abs(theta)<MIN_NUMBER)
							{
							/////////////////////////////////////equation R - theta=0//////////////////////////////////////
							for (unsigned int k=0;k<N;k++)
								{
								if (abs(omega_1(x,k))>MIN_NUMBER)
									{
									unsigned int m=k*NT;
									DDS.coeffRef(2*j,2*m+1) += -2.0*omega_1(x,k);
									minusDS(2*j) 			+= 2.0*omega_1(x,k)*p(2*m+1);
									}
								}
							////////////////////////////////////////////////////////////////////////////////////////
							}
						else
							{
							for (unsigned int k=0;k<N;k++)
								{
								if (abs(omega_1(x,k))>MIN_NUMBER)
									{
									/////////////////////equation I - theta!=0//////////////
									unsigned int m=k*NT;
									DDS.coeffRef(2*j+1,2*m) += (1.0-Gamma)*omega_1(x,k)/(1.0+Gamma);
									minusDS(2*j+1) 			+= -(1.0-Gamma)*omega_1(x,k)*(p(2*m)-minima[0])/(1.0+Gamma);
									/////////////////////equation R - theta!=0//////////////
									minusDS(2*j) 			+= p(2*m+1)*omega_1(x,k)*(1+Gamma)*theta/(1-Gamma);
									DDS.coeffRef(2*j,2*m+1)	+= -omega_1(x,k)*(1.0+Gamma)*theta/(1.0-Gamma);
									bound 				+= -(1.0-Gamma)*omega_1(x,k)*(p(2*j)-minima[0])*(p(2*m)-minima[0])/(1.0+Gamma)\
																 + (1.0+Gamma)*omega_1(x,k)*p(2*j+1)*p(2*m+1)/(1.0-Gamma);
									}
								}
							//////////////////////////////////////equation R - theta!=0//////////////////////////////
				            minusDS(2*j) 					+= theta*Dx*real(Cp(j+1)/dt);
			            	DDS.coeffRef(2*j,2*(j+1)) 		+= -theta*real(Dx/dt);
			            	DDS.coeffRef(2*j,2*(j+1)+1) 	+= theta*imag(Dx/dt);
						    minusDS(2*j) 					+= theta*real(-Dx*Cp(j)/dt);
					    	DDS.coeffRef(2*j,2*j) 			+= theta*real(Dx/dt);
					    	DDS.coeffRef(2*j,2*j+1) 		+= theta*imag(-Dx/dt);
							}
						}
					else ///////////////////////////////// including other terms in action at t=0 ///////////////////////////
						{
						kineticS 	+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						potV 		+= Dt*Dx*V(Cp(j));
						pot_r 		+= Dt*Dx*Vr(Cp(j));
						derivErg(t) += pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						potErg(t) 	+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
						erg(t) 		+= pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0 + Dx*V(Cp(j)) + Dx*Vr(Cp(j));
						//////////////////////////////////////equation I - both///////////////////////////////////
						for (unsigned int k=1; k<2*2; k++)
					    	{
					        int sign = pow(-1,k+1);
					        int direc = (int)(k/2.0);
					        int neighb = neigh(j,direc,sign,NT,N);
					        double dxd = (sign==1? dx: dxm);
					        if (direc == 0)
					        	{
					            minusDS(2*j+1) 					+= Dx*imag(Cp(j+sign)/dt);
					            DDS.coeffRef(2*j+1,2*(j+sign)) 	+= -imag(Dx/dt);
					            DDS.coeffRef(2*j+1,2*(j+sign)+1)+= -real(Dx/dt);
					            }
					        else if (neighb!=-1)
					        	{
					            minusDS(2*j+1) 					+= - imag(Dt*Cp(neighb))/dxd;
					            DDS.coeffRef(2*j+1,2*neighb) 	+= imag(Dt)/dxd;
					            DDS.coeffRef(2*j+1,2*neighb+1) 	+= real(Dt)/dxd;
					            }
					        }
					    comp temp0 = Dx/dt - Dt*(1.0/dx+1.0/dxm);
					    comp temp1 = Dt*Dx*( dV(Cp(j)) + dVr(Cp(j)) );//dV terms should be small
					    comp temp2 = Dt*Dx*(ddV(Cp(j)) + ddVr(Cp(j)));
					    
					    minusDS(2*j+1) 				+= imag(-temp0*Cp(j) + temp1 );
					    DDS.coeffRef(2*j+1,2*j) 	+= imag(temp0 - temp2 );
					    DDS.coeffRef(2*j+1,2*j+1) 	+= real(temp0 - temp2 );
						/////////////////////////////////////////////////////////////////////////////////////////
						if (abs(theta)<MIN_NUMBER)
							{
							//simplest boundary conditions replaced by ones continuously connected to theta!=0 ones
							//DDS.insert(2*j+1,2*(j+1)+1) = 1.0; //zero imaginary part of time derivative
							//DDS.insert(2*j,2*j+1) = 1.0; //zero imaginary part
						
							/////////////////////////////////////equation R - theta=0//////////////////////////////////////
							for (unsigned int k=0;k<N;k++)
								{
								if (abs(omega_1(x,k))>MIN_NUMBER)
									{
									unsigned int m=k*NT;
									DDS.coeffRef(2*j,2*m+1) += -2.0*omega_1(x,k);
									minusDS(2*j) 			+= 2.0*omega_1(x,k)*p(2*m+1);
									}
								}
							////////////////////////////////////////////////////////////////////////////////////////
							}
						else
							{
							for (unsigned int k=0;k<N;k++)
								{
								if (abs(omega_1(x,k))>MIN_NUMBER)
									{
									/////////////////////equation I - theta!=0//////////////
									unsigned int m=k*NT;
									DDS.coeffRef(2*j+1,2*m) += (1.0-Gamma)*omega_1(x,k)/(1.0+Gamma);
									minusDS(2*j+1) 			+= -(1.0-Gamma)*omega_1(x,k)*(p(2*m)-minima[0])/(1.0+Gamma);
									/////////////////////equation R - theta!=0//////////////
									minusDS(2*j) 			+= p(2*m+1)*omega_1(x,k)*(1+Gamma)*theta/(1-Gamma);
									DDS.coeffRef(2*j,2*m+1)	+= -omega_1(x,k)*(1.0+Gamma)*theta/(1.0-Gamma);
									bound 				+= -(1.0-Gamma)*omega_1(x,k)*(p(2*j)-minima[0])*(p(2*m)-minima[0])/(1.0+Gamma)\
																 + (1.0+Gamma)*omega_1(x,k)*p(2*j+1)*p(2*m+1)/(1.0-Gamma);
									}
								}
							//////////////////////////////////////equation R - theta!=0//////////////////////////////
							for (unsigned int k=1; k<2*2; k++)
						    	{
						        int sign = pow(-1,k+1);
						        int direc = (int)(k/2.0);
						        double dxd = (sign==1? dx: dxm);
						        if (direc == 0)
						        	{
						            minusDS(2*j) 					+= Dx*real(Cp(j+sign)/dt)*theta;
					            	DDS.coeffRef(2*j,2*(j+sign)) 	+= -real(Dx/dt)*theta;
					            	DDS.coeffRef(2*j,2*(j+sign)+1) 	+= imag(Dx/dt)*theta;
						            }
						        else
						        	{
						            unsigned int neighb = neigh(j,direc,sign,NT,N);
						            minusDS(2*j) 					+= - real(Dt*Cp(neighb))*theta/dxd;
						            DDS.coeffRef(2*j,2*neighb) 		+= real(Dt)*theta/dxd;
					            	DDS.coeffRef(2*j,2*neighb+1) 	+= -imag(Dt)*theta/dxd;
						            }
						        }
						    minusDS(2*j) 			+= theta*real(-temp0*Cp(j) + temp1 );
					    	DDS.coeffRef(2*j,2*j) 	+= theta*real(temp0 - temp2 );
					    	DDS.coeffRef(2*j,2*j+1) += theta*imag(-temp0 + temp2 );
							}
						}
					}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//bulk
				else
					{
					kineticS 	+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
					kineticT 	+= Dx*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					potV 		+= Dt*Dx*V(Cp(j));
					pot_r 		+= Dt*Dx*Vr(Cp(j));
					erg(t) 		+= Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0 + Dx*V(Cp(j)) + Dx*Vr(Cp(j));
					derivErg(t) += Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
					potErg(t) 	+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
				
		            for (unsigned int k=0; k<2*2; k++)
                	{
                    int sign = pow(-1,k);
                    int direc = (int)(k/2.0);
                    comp dtd = (sign==1? dt: dtm);
                    double dxd = (sign==1? dx: dxm);
                    if (direc == 0)
                    	{
                        minusDS(2*j) 					+= real(Dx*Cp(j+sign)/dtd);
                        minusDS(2*j+1) 					+= imag(Dx*Cp(j+sign)/dtd);
                        DDS.insert(2*j,2*(j+sign)) 		= -real(Dx/dtd);
                        DDS.insert(2*j,2*(j+sign)+1) 	= imag(Dx/dtd);
                        DDS.insert(2*j+1,2*(j+sign)) 	= -imag(Dx/dtd);
                        DDS.insert(2*j+1,2*(j+sign)+1) 	= -real(Dx/dtd);
                        }
                    else
                    	{
                        unsigned int neighb = neigh(j,direc,sign,NT,N);
                        
                        minusDS(2*j) 					+= -real(Dt*Cp(neighb)/dxd);
                        minusDS(2*j+1) 					+= -imag(Dt*Cp(neighb)/dxd);
                        DDS.insert(2*j,2*neighb) 		= real(Dt/dxd);
                        DDS.insert(2*j,2*neighb+1) 		= -imag(Dt/dxd);
                        DDS.insert(2*j+1,2*neighb) 		= imag(Dt/dxd);
                        DDS.insert(2*j+1,2*neighb+1) 	= real(Dt/dxd);
                        }
                    }
		            comp temp0 = Dx*(1.0/dt + 1.0/dtm) - Dt*(1.0/dx + 1.0/dxm);
		            comp temp1 = Dt*Dx*(dV(Cp(j)) + dVr(Cp(j)));
		            comp temp2 = Dt*Dx*(ddV(Cp(j)) + ddVr(Cp(j)));
		                
		            minusDS(2*j) 			+= real(temp1 - temp0*Cp(j));
		            minusDS(2*j+1) 			+= imag(temp1 - temp0*Cp(j));
		            DDS.insert(2*j,2*j) 	= real(-temp2 + temp0);
		            DDS.insert(2*j,2*j+1) 	= imag(temp2 - temp0);
		            DDS.insert(2*j+1,2*j) 	= imag(-temp2 + temp0);
		            DDS.insert(2*j+1,2*j+1) = real(-temp2 + temp0);
		            }
		        }
		    if (pot[0]=='3') DDS.insert(2*N*NT,2*N*NT) = 1.0;
		    action = kineticT - kineticS - potV - pot_r;
		    linErgOffShell(NT-1) = linErgOffShell(NT-2);
	    	linNumOffShell(NT-1) = linNumOffShell(NT-2);
	    	
		    if (pot[0]=='3') {
		    	action 		*= 4.0*pi;
		    	derivErg 	*= 4.0*pi;
		    	potErg 		*= 4.0*pi;
		    	erg			*= 4.0*pi;
		    	linErg		*= 4.0*pi;
		    	linNum		*= 4.0*pi;
		    	linNumOffShell *= 4.0*pi;
		    	linErgOffShell *= 4.0*pi;
		    }
		   
	    	
		    if (abs(theta)<MIN_NUMBER) {
		    	linErg = linErgOffShell;
		    	linNum = linNumOffShell;
		    }
		    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			//checking linErg, linNum, bound, computing W, E
			
			//trivial test that erg=potErg+derivErg
			if (false) {
				for (unsigned int j=0; j<NT; j++) {
					double diff = absDiff(erg(j),potErg(j)+derivErg(j));
					if (diff>1.0e-14) cerr << "erg(" << j << ") != potErg + derivErg. absDiff = " << diff << endl;
				}
			}
			
			//calculating continuum approx to linErg and linNum on initial time slice
			if (pot[0]=='3') {
				for (unsigned int k=1; k<cutoff; k++) {
					double momtm = k*pi/(L-r0);
					double freqSqrd = 1.0+pow(momtm,2.0);
					double Asqrd, integral1 = 0.0, integral2 = 0.0;
					for (unsigned int l=0; l<N; l++) {
						double r = r0 + l*a;
						unsigned int m = l*NT;
						integral1 += a*p(2*m)*pow(2.0/L,0.5)*sin(momtm*r);
						integral2 += a*(p(2*(m+1))-p(2*m))*pow(2.0/L,0.5)*sin(momtm*r)/b;
					}
					Asqrd = pow(integral1,2.0) + pow(integral2,2.0)/freqSqrd;
					linErgContm += 2.0*pi*Asqrd*freqSqrd;
					linNumContm += 2.0*pi*Asqrd*pow(freqSqrd,0.5);
				}
			}
			
			//calculating a_k and b*_k
			cVec a_k(N), b_k(N); //N.b. b_k means b*_k but can't put * in the name
			a_k = Eigen::VectorXcd::Zero(N);
			b_k = Eigen::VectorXcd::Zero(N);
			double T0 = real(coord(0,0));
			comp dt0 = dtFn(0);
			for (unsigned int n=0; n<N; n++) {
				double w_n = sqrt(eigenValues(n));
				for (unsigned int j=0; j<N; j++) {
					unsigned int m=j*NT;
					double sqrtDj = sqrt(4.0*pi*DxFn(j));
					if (abs(eigenValues(n))>1.0e-16) {
						a_k(n) += exp(ii*w_n*T0)*sqrt(2.0*w_n)*eigenVectors(j,n)* \
									sqrtDj*((Cp(m+1)-minima[0])-(Cp(m)-minima[0])*exp(ii*w_n*dt0)) \
										/(exp(-ii*w_n*dt0)-exp(ii*w_n*dt0));
						b_k(n) += exp(-ii*w_n*T0)*sqrt(2.0*w_n)*eigenVectors(j,n)* \
									sqrtDj*((Cp(m+1)-minima[0])-(Cp(m)-minima[0])*exp(-ii*w_n*dt0)) \
										/(exp(ii*w_n*dt0)-exp(-ii*w_n*dt0));
					}
				}
			}			
			
			//using a_k and b*_k to check that a_k=Gamma*b*_k and p->p_lin as t->0 and that Sum_k(a_k*b*_k)=linNum
			// and that Sum_k(w_k*a_k*b*_k)=linErg
			double ABtest = 0.0, linRepTest, ABNEtest;
			linNumAB = 0.0, linErgAB = 0.0;
			cVec linRep(N), p0(N);
			linRep = Eigen::VectorXcd::Constant(N,minima[0]);
			for (unsigned int j=0; j<N; j++) {
				ABtest += absDiff(a_k(j),conj(b_k(j))*Gamma);
				p0(j) = Cp(j*NT);
				if (j<cutoff) {
					linNumAB += a_k(j)*b_k(j);
					linErgAB += sqrt(eigenValues(j))*a_k(j)*b_k(j);
				}
				for (unsigned int n=0; n<N; n++) {
					double w_n = sqrt(eigenValues(n));
					double sqrtDj = sqrt(4.0*pi*DxFn(j));
					if (abs(w_n)>1.0e-16) {
						linRep(j) += eigenVectors(j,n)*(a_k(n)*exp(-ii*w_n*T0)+b_k(n)*exp(ii*w_n*T0)) \
										/sqrt(2.0*w_n)/sqrtDj;
					}
				}
			}
			ABtest /= (double)N;
			linRepTest = absDiff(p0,linRep);
			double ABNtest = absDiff(linNumAB,linNum(0));
			double ABEtest = absDiff(linErgAB,linErg(0));
			ABNEtest = (ABNtest>ABEtest? ABNtest: ABEtest);
			AB_test.push_back(ABtest);
			linRep_test.push_back(linRepTest);
			ABNE_test.push_back(ABNEtest);
			/*if (AB_test.back()>closenessAB || !isfinite(AB_test.back())) cerr << "AB_test = " << AB_test.back() << endl;
			if (linRep_test.back()>closenessLR || !isfinite(linRep_test.back()))cerr << "linRep_test = " << linRep_test.back() << endl;
			if (ABNE_test.back()>closenessABNE || !isfinite(ABNE_test.back())) {
				cerr << "linNumAB = " << linNumAB << ", linNum(0) = " << \
					linNum(0) << ", linNumOffShell(0) = " << linNumOffShell(0) << endl;
				cerr << "linErg = " << linErgAB << ", linErg(0) = " << linErg(0) \
					<< ", linErgOffShell(0) = " << linErgOffShell(0) << endl;
			}*/
			
			//checking pot_r is much smaller than the other potential terms
			reg_test.push_back(abs(pot_r/potV));
			if (reg_test.back()>closenessR)
				{
				cerr << "regularisation term is too large, regTest = " << reg_test.back() << endl;
				}
			
			//checking linearisation of linErg and linNum
			double linTestE;	double linTestN;
			double linEMax = 0.0;	double linEMin = 5.0e15; //surely it's going to be less than this
			double linNMax = 0.0;	double linNMin = 5.0e15;
			unsigned int linearInt = (unsigned int)(Na/10);
			for (unsigned int j=1;j<(linearInt+1);j++)
				{
				if (abs(linErgOffShell(j))>linEMax) linEMax = abs(linErgOffShell(j));
				if (abs(linErgOffShell(j))<linEMin) linEMin = abs(linErgOffShell(j));
				if (abs(linNumOffShell(j))>linNMax) linNMax = abs(linNumOffShell(j));
				if (abs(linNumOffShell(j))<linNMin) linNMin = abs(linNumOffShell(j));
				}
			linTestE = absDiff(linEMax,linEMin);
			linTestN = absDiff(linNMax,linNMin);
			lin_test.push_back(linTestE);
			if (linTestN>linTestE) lin_test.back() = linTestN;
			
			//checking agreement of on and off shell linear energy at initial time
			double onShellTest = absDiff(linErg(0),linErgOffShell(0));
			onShell_test.push_back(onShellTest);
			
			//checking conservation of E
			double conservTest = absDiff(erg(1),erg(NT-2));
			conserv_test.push_back(conservTest);
						
			//defining E, Num and cW
			E = real(linErg(0));
			Num = real(linNum(0));
			E_exact = 0.0;
			for (unsigned int j=1; j<(linearInt+1); j++) E_exact += real(erg(j));
			E_exact /= (double)linearInt;
			W = - E*2.0*Tb - theta*Num - bound + 2.0*imag(action);
			
			//testing imaginary part of energy
			double imErgTest = 0.0;
			for (unsigned int j=0; j<NT; j++) if (abs(erg(j))>MIN_NUMBER) imErgTest += imag(erg(j))/abs(erg(j));
			imErgTest /= (double)NT;
			imErg_test.push_back(imErgTest);
			
			//checking agreement between erg and linErg
			double trueTest = absDiff(E,E_exact);
			true_test.push_back(trueTest);
			if (!isfinite(trueTest)) cout << "E = " << E << ", E_exact = " << E_exact << ", linearInt = " << linearInt << endl;
			
			//checking lattice small enough for E, should have parameter for this
			double momTest = E*b/Num/pi; //perhaps should have a not b here
			mom_test.push_back(momTest);
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			//solving for delta in DDS*delta=minusDS, where p' = p + delta		
			vec delta(2*N*NT+2);
			delta = Eigen::VectorXd::Zero(2*N*NT+2);
			DDS.prune(MIN_NUMBER);
			DDS.makeCompressed();
			Eigen::SparseLU<spMat> solver;
		
			solver.analyzePattern(DDS);
			if(solver.info()!=Eigen::Success)
				{
				cerr << "DDS pattern analysis failed, solver.info() = "<< solver.info() << endl;
				return 1;
				}		
			solver.factorize(DDS);
			if(solver.info()!=Eigen::Success) 
				{
				cerr << "Factorization failed, solver.info() = "<< solver.info() << endl;
				return 1;
				}
			delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side
			if(solver.info()!=Eigen::Success)
				{
				cerr << "Solving failed, solver.info() = "<< solver.info() << endl;
				cerr << "log(abs(det(DDS))) = " << solver.logAbsDeterminant() << endl;
				cerr << "sign(det(DDS)) = " << solver.signDeterminant() << endl;
				return 1;
				}
		
			//independent check on whether calculation worked
			vec diff(2*N*NT+2);
			diff = DDS*delta-minusDS;
			double maxDiff = diff.maxCoeff();
			maxDiff = abs(maxDiff);
			calc_test.push_back(maxDiff);
			if (calc_test.back()>closenessC)
				{
				cerr << "Calculation failed" << endl;
				cerr << "calc_test = " << calc_test.back() << endl;
				return 1;
				}

			//assigning values to phi
			p += delta;
		
			//passing changes on to complex vector
			Cp = vecComplex(p,N*NT);
			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			string prefix = "./data/" + timeNumber;
			string suffix = "_" + numberToString<unsigned int>(fileLoop) + "_" + numberToString<unsigned int>(loop)+"_"+numberToString<unsigned int>(runs_count)+".dat";
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
			if ((print_choice.compare("d")==0 || print_choice.compare("e")==0))
				{
				string earlydeltaFile = prefix + "maindeltaE" + suffix;
				printVector(earlydeltaFile,delta);
				}
			if ((print_choice.compare("l")==0 || print_choice.compare("e")==0))
				{
				string earlyLinErgFile = prefix + "mainlinErgE"+suffix;
				simplePrintCVector(earlyLinErgFile,linErg);
				//gpSimple(earlyLinErgFile);
				string earlyErgFile = prefix + "mainergE" + suffix;
				simplePrintCVector(earlyErgFile,erg);
				//gpSimple(earlyErgFile);
				string earlyDerivErgFile = prefix + "mainderivErgE"+suffix;
				//derivErg.conservativeResize(Na);
				simplePrintCVector(earlyDerivErgFile,derivErg);
				string earlyPotErgFile = prefix + "mainpotErgE"+suffix;
				//potErg.conservativeResize(Na);
				simplePrintCVector(earlyPotErgFile,potErg);
				}
			if ((print_choice.compare("ab")==0 || print_choice.compare("e")==0))
				{
				string earlya_kFile = prefix + "maina_kE"+suffix;
				simplePrintCVector(earlya_kFile,a_k);
				string earlyb_kFile = prefix + "mainb_kE" + suffix;
				simplePrintCVector(earlyb_kFile,b_k);
				}
			}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			//convergence issues
			//evaluating norms
			double normDS = minusDS.norm();
			double maxDS = minusDS.maxCoeff();
			double minDS = minusDS.minCoeff();
			if (-minDS>maxDS) maxDS = -minDS;
			double maxP = p.maxCoeff();
			double minP = p.minCoeff();
			if (-minP>maxP) maxP = -minP;
			double normP = p.norm();
			double normDelta = delta.norm();
		
			//assigning test values
			//quantities used to stop newton-raphson loop
			action_test.push_back(absDiff(action,action_last));
			action_last = action;
			sol_test.push_back(normDS/normP);
			solM_test.push_back(maxDS/maxP);
			delta_test.push_back(normDelta/normP);
			
			//printing tests to see convergence
			if (runs_count==1)
				{
				printf("%5s%5s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s\n","loop","run","solTest","solMTest","deltaTest","linTest","trueTest","shellTest","ABTest","ABNETest","cnsrvTest","momTest");
				}
			printf("%5i%5i%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g\n",loop,runs_count,sol_test.back(),solM_test.back(),delta_test.back(),lin_test.back(),true_test.back(),onShell_test.back(),AB_test.back(),ABNE_test.back(),conserv_test.back(),mom_test.back());
			
			if (delta_test.back()>closenessD) {
				cerr << "deltaTest > " << closenessD << ", stopping NR loop." << endl;
				break;
			}
			
			} //ending while loop
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	    		//misc end of program tasks - mostly printing

		//checking different measures of energy agree
		if (true_test.back()>closenessT || onShell_test.back()>closenessON || ABNE_test.back()>closenessABNE) {
			cout << "erg(1)             = " << erg(1) << endl;
			cout << "linErgOffShell(0)  = " << linErgOffShell(0) << endl;
			cout << "linErgAB           = " << linErgAB << endl;
			cout << "linErgContm        = " << linErgContm << endl;
			cout << "linErg(0)          = " << linErg(0) << endl;
		}

		//checking energy conserved
		if (conserv_test.back()>closenessCon) cout << endl << "ergTest = " << conserv_test.back() << endl;
		
		//checking energy real
		if (imErg_test.back()>closenessIE) cout << endl << "imErg = "<< imErg_test.back()  << endl;
		
		//checking lattice small enough
		if (mom_test.back()>closenessP) cout << endl << "momTest = "<< mom_test.back()  << endl;
		
		if (pot[0]=='3' && abs(theta)<MIN_NUMBER) {
			comp temp = linErg(0);
			if (absDiff(temp,linErgContm)>closenessCL) 
				cout << "linErg(0) = " << temp << " and linErgContm = " << linErgContm << " don't agree. E_exact = " << E_exact << endl;
			temp = linNum(0);
			if (absDiff(temp,linNumContm)>closenessCL) 
				cout << "linNum(0) = " << temp << " and linNumContm = " << linNumContm << " don't agree. " << endl;
		}
		
		//stopping clock
		time = clock() - time;
		double realtime = time/1000000.0;
	
		//printing to terminal
		printf("\n");
		printf("%8s%8s%8s%8s%8s%8s%8s%8s%14s%14s%14s%14s\n","runs","time","N","NT","L","Tb","dE","theta","Num","E","im(action)","W");
		printf("%8i%8g%8i%8i%8g%8g%8g%8g%14g%14g%14g%14g\n",runs_count,realtime,N,NT,L,Tb,dE,theta,Num,E,imag(action),real(W));
		printf("\n");
		 printf("%60s\n","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

		//printing action value
		FILE * actionfile;
		actionfile = fopen("./data/mainAction.dat","a");
		fprintf(actionfile,"%12s%6i%6i%8g%8g%8g%8g%14g%14g%12g%14g%14g%14g\n",timeNumber.c_str(),N,NT,L,Tb,dE,theta,E,Num\
		,real(W),sol_test.back(),lin_test.back(),true_test.back());
		fclose(actionfile);
		
		string prefix = "./data/" + timeNumber;
		string suffix = "_" + numberToString<unsigned int>(fileLoop)+"_" + numberToString<unsigned int>(loop)+ ".dat";
	
		//copying a version of inputs with timeNumber and theta changed
		string runInputs = prefix + "inputsM_"+ numberToString<unsigned int>(fileLoop) + "_" + numberToString<unsigned int>(loop); //different suffix
		if (abs(maxTheta-minTheta)>MIN_NUMBER) 	changeInputs(runInputs,"theta", numberToString<double>(theta),inputsFiles[fileLoop]);
		else if (abs(maxTb-minTb)>MIN_NUMBER)	changeInputs(runInputs,"Tb", numberToString<double>(Tb),inputsFiles[fileLoop]);
		else 									copyFile(inputsFiles[fileLoop],runInputs);
		printf("%12s%30s\n","output: ",runInputs.c_str());
	
		//printing output phi
		string tpifile =  prefix + "mainpi"+suffix;
		printVector(tpifile,p);
		//gp(tpifile,"repi.gp");
		printf("%12s%30s\n"," ",tpifile.c_str());
	
		//printing output minusDS				
		string minusDSfile = prefix + "mainminusDS"+suffix;
		printVector(minusDSfile,minusDS);
		printf("%12s%30s\n"," ",minusDSfile.c_str());
				
		//printing output DDS
		string DDSfile = prefix + "mainDDS"+suffix;
		printSpmat(DDSfile,DDS);
	    printf("%12s%30s\n"," ",DDSfile.c_str());
	    
		//printing linNum
		string linNumFile = prefix + "mainlinNum"+suffix;
		linNum.conservativeResize(Na);
		simplePrintCVector(linNumFile,linNum);
		printf("%12s%30s\n"," ",linNumFile.c_str());
		//gpSimple(linNumFile);
	
		//printing linErg
		string linErgFile = prefix + "mainlinErg"+suffix;
		linErg.conservativeResize(Na);
		simplePrintCVector(linErgFile,linErg);
		//gpSimple(linErgFile);
		printf("%12s%30s\n"," ",linErgFile.c_str());
		
		//printing linErgOffShell
		string linErgOffShellFile = prefix + "mainlinErgOffShell"+suffix;
		linErgOffShell.conservativeResize(Na);
		simplePrintCVector(linErgOffShellFile,linErgOffShell);
		//gpSimple(linErgFileOffShell);
		printf("%12s%30s\n"," ",linErgOffShellFile.c_str());
		
		//printing derivErg
		string derivErgFile = prefix + "mainderivErg"+suffix;
		//derivErg.conservativeResize(Na);
		simplePrintCVector(derivErgFile,derivErg);
		//gpSimple(linErgFile);
		printf("%12s%30s\n"," ",derivErgFile.c_str());
		
		//printing potErg
		string potErgFile = prefix + "mainpotErg"+suffix;
		//potErg.conservativeResize(Na);
		simplePrintCVector(potErgFile,potErg);
		//gpSimple(linErgFile);
		printf("%12s%30s\n"," ",potErgFile.c_str());
	
		//printing erg
		string ergFile = prefix + "mainerg" + suffix;
		erg.conservativeResize(Na);
		simplePrintCVector(ergFile,erg);
		//gpSimple(ergFile);
		printf("%12s%30s\n"," ",ergFile.c_str());
		cout << endl;
		
		} //ending parameter loop
	} //ending file loop

return 0;
}

