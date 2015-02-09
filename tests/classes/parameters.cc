/*
	definitions for functions and classes for dealing with parameters
*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
parameter related errors
	- unset
	- load
-------------------------------------------------------------------------------------------------------------------------*/

// unset
string ParameterError::Unset::message() const{
	return "Parameter " + pName + " not set.";
}

// load
string ParameterError::Load::message() const{
	return "Error: parameters not loaded from " + filename + ".";
}

/*-------------------------------------------------------------------------------------------------------------------------
PrimaryParameters member functions
	- save
	- load
-------------------------------------------------------------------------------------------------------------------------*/

//save
void PrimaryParameters::save(const string& filename) const {
	ofstream os;
	os.open(filename.c_str());
	if (!os.good()) {
		FileError::StreamNotGood e(filename);
		throw e;
	}
	os << left;
	os << setw(20) << "pot" << setw(20) << pot << endl;
	os << setw(20) << "N" << setw(20) << N << endl;
	os << setw(20) << "Na" << setw(20) << Na << endl;
	os << setw(20) << "Nb" << setw(20) << Nb << endl;
	os << setw(20) << "Nc" << setw(20) << Nc << endl;
	os << setw(20) << "LoR" << setw(20) << LoR << endl;
	os << setw(20) << "dE" << setw(20) << dE << endl;
	os << setw(20) << "Tb" << setw(20) << Tb << endl;
	os << setw(20) << "theta" << setw(20) << theta << endl;
	os << setw(20) << "reg" << setw(20) << reg << endl;
	os << endl;
	os.close();
}

//load
void PrimaryParameters::load(const string& filename) {
	ifstream is;
	is.open(filename.c_str());
	if (!is.good()) {
		FileError::StreamNotGood e(filename);
		throw e;
	}
	string line;
	bool set = false;
	while(getline(is,line)) {
		if (line[0]=="#") continue;
		if (line.empty()) break;
		set = true;
		ss >> "pot" >> pot;
		ss >> "N" >> N;
		ss >> "Na" >> Na;
		ss >> "Nb" >> Nb;
		ss >> "Nc" >> Nc;
		ss >> "LoR" >> LoR;
		ss >> "dE" >> dE;
		ss >> "Tb" >> Tb;
		ss >> "theta" >> theta;
		ss >> "reg" >> reg;
		ss.close();
	}
	if (!set) {
		ParameterError::Load e(filename);
		throw e;
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
SecondaryParameters member functions
	- set SecondaryParameters from PrimaryParameters
-------------------------------------------------------------------------------------------------------------------------*/

// set secondary parameters
void SecondaryParameters::setSecondaryParameters (const struct PrimaryParameters&) {
	NT = Na + Nb + Nc;
	
	double epsilon;						// determined by dE
	double R; 							// size of bubble
	double Gamma; 						// equals exp(-theta)
	double L;
	double a; 							// step sizes in each spatial dimension
	double b; 							// step sizes in time
	double Ta;
	double Tc;
	vector<double> minima(2);			// minima of V
	double mass2; 						// as derived from V''
}

/*-------------------------------------------------------------------------------------------------------------------------
Parameters member functions
	- empty constructor
	- constructor from PrimaryParameters and SecondaryParameters
	- print to shell
	- update secondary parameters
	- change parameters based on change in one
-------------------------------------------------------------------------------------------------------------------------*/

Parameters::Parameters(): PrimaryParmaeters(), SecondaryParameters()	{}			// empty constructor

Parameters::Parameters(const PrimaryParameters& p1, const SecondaryParmeters& p2): \
					PrimaryParmaeters(p1), SecondaryParameters(p2)	{}				// constructor using primary and secondary parameters

void Parameters::printParameters() const {
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","N","Na","Nb","Nc","L","Ta","Tb","R","dE","theta","reg");
	printf("%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g%8g\n",\
			N,Na,Nb,Nc,\
			L,Ta,Tb,R,\
			dE,theta,reg);
	printf("\n");
}
