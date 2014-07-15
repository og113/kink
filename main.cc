#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "pf.h"

//struct to hold answers to questions
struct aq{
	unsigned int fileNo;
	double maxTheta;
	unsigned int totalLoops;
};

int main()
{
//asking questions
aq Xaq;
cout << "'which data/pic#.dat file to load? (#) ";
cin >> Xaq.fileNo;
cout << endl << "input final value of angle (input 0 for no steps) ";
cin >> Xaq.maxTheta;
cout << endl << "input number of loops to get there";
cin >> Xaq.totalLoops;

//loading initial phi
string inFile = "./data/pic"+to_string(Xaq.fileNo)+".dat";
ofstream f;
f.open((inFile).c_str());
///////////
f.close();

