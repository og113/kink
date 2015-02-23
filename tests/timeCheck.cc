//program to check time coord
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include "gnuplot_i.hpp"
#include "pf.h"
#include "files.h"

using namespace std;

int main()
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//getting variables and user inputs from inputs

//defining the time to label output
string timeNumber = currentDateTime();

ifstream fin;
fin.open("inputs");
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
			if (absolute(theta)>2.0e-16)
				{
				cout << "theta != 0" << endl;
				cout << "theta = " << theta << endl;
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
	cout << "unable to open inputs" << endl;
	}
fin.close();
inP = aq.inputChoice; //just because I write this a lot
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//calculated quantities derived from the inputs, and loading eigVec and eigVal	

//derived quantities
NT = Na + Nb + Nc;
Gamma = exp(-theta);
epsilon = dE;
R = 2.0/3.0/epsilon;
alpha *= R;
L = LoR*R;
if (inP.compare("p") == 0 || inP.compare("f") == 0)
	{
	if (Tb<R)
		{
		angle = asin(Tb/R);
		double Ltemp = 1.5*(1.5*Tb*tan(angle));
		if (Ltemp<L) //making sure to use the smaller of the two possible Ls
			{
			L=Ltemp;
			}
		}
	}
else if (inP.compare("b") == 0)
	{
	Tb = 1.5*R;
	}
a = L/(N-1.0);
b = Tb/(Nb-1.0);
Ta = b*Na;
Tc = b*Nc;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//defining time coord t
cVec t(NT);
for (unsigned int j=0;j<NT;j++)
	{
	t(j) = simpleTime(j);
	}
string filename = "./data/timeCheck.dat";
simplePrintCVector(filename,t);

gpSimple2(filename);

return 0;
}
