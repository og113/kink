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
	string printChoice;
	unsigned int printRun;
};

int main()
{
//asking questions
aq Xaq;
cout << "'which data/picOut#.mat file to load? (#) ";
cin >> Xaq.

