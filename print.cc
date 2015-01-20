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

int main(int argc, char** argv)
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//putting together filename from questions

string filename, id, method;
if (argc==1) {
	cout << "filename/n: ";
	cin >> filename;
}
else if (argc==2) {
	filename = argv[1];
}
else if (argc==3) {
	filename = argv[1];
	string temp = argv[2];
	if (temp[0]=='-') temp = temp.substr(1);
	method = temp;
}
if (filename.compare("n")==0)
	{
	filename = "";
	string timeNumber;
	cout << "timeNumber: ";
	cin >> timeNumber;
	char program;
	cout << "pi or main, p/m: ";
	cin >> program;
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
	if (program=='m')
		{
		filename  = "./data/" + timeNumber + "main" + id + "_" + to_string(fileLoop) + "_" + to_string(paramLoop);
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
	}
else if (method.empty())
	{
	cout << "printing method: simple or vector (s/v): ";
	cin >> method;
	cout << endl; 
	}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//loading parameters

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//printing

if (id.compare("linErg")==0 || id.compare("linNum")==0 || id.compare("erg")==0 || method.compare("s")==0)
	{
	gpSimple(filename);
	}
else if (id.compare("tpi")==0 || id.compare("pi")==0 || id.compare("chiX")==0 || id.compare("chiT")==0 || method.compare("v")==0)
	{
	gp(filename,"repi.gp");
	}
else if (id.compare("minusDS")==0)
	{
	gp(filename,"pi.gp");
	}
else
	{
	cout << "print does not recognise input" << endl;
	}

return 0;
}
