//program to speed up printing via gnuplot
#include <Eigen/Sparse>
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
#include <gsl/gsl_poly.h>
#include "gnuplot_i.hpp"
#include "pf.h"

using namespace std;

int main()
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//putting together filename from questions

string timeNumber;
cout << "timeNumber: ";
cin >> timeNumber;
char program;
cout << "pi or main, p/m: ";
cin >> program;
string id;
cout << "identifier: ";
cin >> id;
unsigned int fileLoop, paramLoop, runsCount;
if (program=='m')
	{
	cout << "fileLoop: ";
	cin >> fileLoop;
	}
cout << "parameter loop: ";
cin >> paramLoop;
string filename;
if (program=='m')
	{
	filename  = "./data/" + timeNumber + "main" + id + to_string(fileLoop) + "_" + to_string(paramLoop);
	}
else if (program=='p')
	{
	filename  = "./data/" + timeNumber + id + "_" + to_string(paramLoop);
	}
else
	{
	cout << "error: program must be either p/m" << endl;
	return 0;
	}
if (id.back()=='p' || id.back()=='b')
	{
	id = id.substr(0,id.size()-1);
	}
if (id.back()=='E')
	{
	cout << "runsCount: ";
	cin >> runsCount;
	filename = filename + "_" + to_string(runsCount);
	id = id.substr(0,id.size()-1);
	}
filename = filename + ".dat";
cout << "using filename: " << filename << endl;
cout << "trimmed id: " << id << endl;
cout << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//loading parameters

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//printing

if (id.compare("linErg")==0 || id.compare("linNum")==0 || id.compare("erg")==0)
	{
	gpSimple(filename);
	}
else if (id.compare("tpi")==0 || id.compare("pi")==0 || id.compare("minusDS")==0|| id.compare("chiX")==0|| id.compare("chiT")==0)
	{
	gp(filename,"repi.gp");
	}
else
	{
	cout << "print does not recognise input" << endl;
	}

return 0;
}
