/*
program to solve the scalar O(N) symmetric bounce equation:
	f'' + (N-1)f'/r - f + f^3 = 0
input argument argv[1]=N
*/
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
#include <ctype.h>
#include <limits> //for DBL_MIN etc
#include <cstring> //for memcpy
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "gnuplot_i.hpp"

using namespace std;

#define pi 3.14159265359

struct paramsStruct {double N;};
struct paramsStruct parameters;

int func (double t, const double y[], double f[], void *params)
	{
	struct paramsStruct * par = (struct paramsStruct *)params;
	double N = (par->N);
	f[0] = y[1];
	f[1] = -(N-1.0)*y[1]/t + y[0] - pow(y[0],3.0);
	f[2] = y[3];
	f[3] = -(N-1.0)*y[3]/t + y[2] - 3.0*y[2]*pow(y[0],2.0);
	return GSL_SUCCESS;
	}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
	{
	struct paramsStruct * par = (struct paramsStruct *)params;
	double N = (par->N);
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 0, 2, 0.0);
	gsl_matrix_set (m, 0, 3, 0.0);
	gsl_matrix_set (m, 1, 0, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 1, 1, -(N-1.0)/t);
	gsl_matrix_set (m, 1, 2, 0.0);
	gsl_matrix_set (m, 1, 3, 0.0);
	gsl_matrix_set (m, 2, 0, 0.0);
	gsl_matrix_set (m, 2, 1, 0.0);
	gsl_matrix_set (m, 2, 2, 0.0);
	gsl_matrix_set (m, 2, 3, 1.0);
	gsl_matrix_set (m, 3, 0, -6.0*y[2]*y[0]);
	gsl_matrix_set (m, 3, 1, 0.0);
	gsl_matrix_set (m, 3, 2, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 3, 3, -(N-1.0)/t);
	dfdt[0] = 0.0;
	dfdt[1] = (N-1.0)*y[1]/pow(t,2.0);
	dfdt[2] = 0.0;
	dfdt[3] = (N-1.0)*y[3]/pow(t,2.0);
	return GSL_SUCCESS;
	}
	
//gives absolute value of a number
double absolute (const double& amplitude)
	{
	double abs_amplitude;
	if (amplitude > 0)
		{
		abs_amplitude = amplitude;
		}
	else
		{
		abs_amplitude = -amplitude;
		}
	return abs_amplitude;
	}

//simple function to print a vector to file	
void simplePrintVector(const string& printFile, const vector<double>& vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F.precision(16);
	F << left;
	unsigned int length = vecToPrint.size();
	for (unsigned int j=0; j<length; j++)
		{
		F << setw(25) << j << setw(25) << vecToPrint[j] << endl;
		}
	F.close();
	}
	
//print repi via gnuplot
void gpSimple(const string & readFile, const string & outFile) 
	{
	string commandOpenStr = "gnuplot -persistent";
	const char * commandOpen = commandOpenStr.c_str();
	FILE * gnuplotPipe = popen (commandOpen,"w");
	string commandStrPic1 = "set term png size 1600,800";
    string commandStrPic2 = "set output '" + outFile + "'";
	string commandStr1 = "plot \"" + readFile + "\" using 1:2 with linespoints";
	string commandStr2 = "pause -1";
	const char * commandPic1 = commandStrPic1.c_str();
	const char * commandPic2 = commandStrPic2.c_str();
	const char * command1 = commandStr1.c_str();
	const char * command2 = commandStr2.c_str();
	fprintf(gnuplotPipe, "%s \n", commandPic1);
	fprintf(gnuplotPipe, "%s \n", commandPic2);
	fprintf(gnuplotPipe, "%s \n", command1);
	fprintf(gnuplotPipe, "%s \n", command2);
	pclose(gnuplotPipe);
	}
	
//getting the date and time
const string currentDateTime()
	{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%y%m%d%H%M%S", &tstruct);
    return buf;
	}
	
int main(int argc,char **argv)
{
if (argc!=2)
	{
	printf("\nargc = %i. Must give argument N.\n\n",argc);
	return 0;
	}
double N = atof(argv[1]);
printf("\nFinding the O(%g) bounce.\n\n",N);
parameters = {N};

gsl_odeiv2_system sys = {func, jac, 4, &parameters};

double F = 1.0, dF;
double aim = 0.0;
double closeness = 1.0e-6;
double r0 = 1.0e-12, r1 = 10.0;
const unsigned int Nx = 1e4;
double dr = r1-r0;
dr /= (double)Nx;
vector<double> y0Vec(Nx+1), y2Vec(Nx+1);
unsigned int runsCount = 0;

double Y0, Y1 = 0.0, Y2 = 1.0, Y3 = 0.0;
if (absolute(N-3.0)<0.1) { Y0 = 4.337;}
else if (absolute(N-2.0)<0.1) { Y0 = 2.206;}
else if (absolute(N-4.0)<0.1) { cout << "initial y0: "; cin >> Y0;}
else { cout << "initial y0: "; cin >> Y0;}

printf("%16s%16s%16s%16s%16s%16s%16s%16s\n","run","Nx","y(r1)","yMin","F-aim","Y0Old","Y0New","-F/dF");

while (absolute(F-aim)>closeness)
	{
	runsCount++;
	
	//if possible, calculating first step, from r=0 to r=dr, using analytic approximation for small r
	bool firstStep = false;
	/*
	if (dr<0.1*pow(absolute(1.0-3*pow(Y0,2.0)),0.5))
		{
		firstStep = true;
		y0Vec[0] = Y0;
		y2Vec[0] = Y2;
		Y1 += dr*(Y0-pow(Y0,3.0))/N;
		Y0 += pow(dr,2.0)*(Y0-pow(Y0,3.0))/2.0/N;
		double temp = Y2;
		Y2 += dr*Y3; //using simple forward euler derivative
		Y3 += dr*(-(N-1.0)*Y3/dr + temp - 3.0*Y3*temp);
		r0 = dr;
		}
	*/
	
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_yp_new (&sys,  gsl_odeiv2_step_rk8pd, 1.0e-8, 1.0e-6, 0.0);

	double r, ri;
	double y0[4] = { Y0, Y1, Y2, Y3};
	double y[4];
	memcpy(y, y0, sizeof(y0));
	double yMin = absolute(y[0]-aim);
	int status;
	unsigned int i, iStart, iMin = 0;
	if (firstStep)
		{
		iStart = 2;
		r = dr;
		y0Vec[1] = y[0];
		y2Vec[1] = y[2];
		}
	else
		{
		iStart = 1;
		r = r0;
		y0Vec[0] = y[0];
		y2Vec[0] = y[2];
		}

	for (i = iStart; i <= Nx; i++)
		{
		ri =(double)i;
		ri *= dr;
		ri += r0;
		status = gsl_odeiv2_driver_apply (d, &r, ri, y);
		if (status != GSL_SUCCESS)
			{
			printf ("error, return value=%d\n", status);
			printf ("i = %3i, r = %3g\n",i,r);
			break;
			}
		if ((y[0]-aim)<yMin && (y[0]-aim)>0.0)
			{
			yMin = y[0]-aim;
			iMin = i;
			}
		y0Vec[i] = y[0];
		y2Vec[i] = y[2];
		//printf ("%.5e %.5e %.5e\n", r, y[0], y[1]);
		if ((y[0]-aim)<0.0)
			{
			iMin = i;
			if ((y[0]-aim)<-0.2) { break; }
			}
		}
	if (status != GSL_SUCCESS){break;}
		
	F = y0Vec[iMin]-aim; //as final boundary condition is y(r1)=0.0;
	dF = y2Vec[iMin];
	printf("%16i%16i%16g%16g%16g%16.12g",runsCount,i,y[0],yMin,F-aim,Y0);
	if (absolute(dF)>2.0e-16)
		{
		if ((Y0-F/dF)>0.0)
			{
			Y0 += -F/dF;
			}
		else
			{
			Y0 /= 1.2; //arbitrary fudge
			}
		printf("%16.12g%16g\n",Y0,-F/dF);
		}
	gsl_odeiv2_driver_free (d);
	if (i==(Nx+1))
		{
		F = y[0];
		}
	}

//printing solution
string filename = "./data/O" + to_string((int)N) + "bounce.dat", picname = "./pics/O" + to_string((int)N) + "bounce.png";
simplePrintVector(filename,y0Vec);
printf("\nSolution printed: %20s  %20s\n\n",filename.c_str(),picname.c_str());
gpSimple(filename,picname);

return 0;
}
