//parameters and functions for pi.cc
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

struct void_params {};
struct void_params paramsVoid;

#define pi 3.14159265359

int func (double t, const double y[], double f[], void *params)
	{
	f[0] = y[1];
	f[1] = -2.0*y[1]/t + y[0] - pow(y[0],3.0);
	f[2] = y[3];
	f[3] = -2.0*y[3]/t + y[2] - 3.0*y[2]*pow(y[0],2.0);
	return GSL_SUCCESS;
	}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
	{
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 0, 2, 0.0);
	gsl_matrix_set (m, 0, 3, 0.0);
	gsl_matrix_set (m, 1, 0, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 1, 1, -2.0/t);
	gsl_matrix_set (m, 1, 2, 0.0);
	gsl_matrix_set (m, 1, 3, 0.0);
	gsl_matrix_set (m, 2, 0, 0.0);
	gsl_matrix_set (m, 2, 1, 0.0);
	gsl_matrix_set (m, 2, 2, 0.0);
	gsl_matrix_set (m, 2, 3, 1.0);
	gsl_matrix_set (m, 3, 0, -6.0*y[2]*y[0]);
	gsl_matrix_set (m, 3, 1, 0.0);
	gsl_matrix_set (m, 3, 2, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 3, 3, -2.0/t);
	dfdt[0] = 0.0;
	dfdt[1] = 2.0*y[1]/pow(t,2.0);
	dfdt[2] = 0.0;
	dfdt[3] = 2.0*y[3]/pow(t,2.0);
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
	
void simplePrintVector(const string& printFile, vector<double> vecToPrint)
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
void gpSimple(const string & readFile) 
	{
	string commandOpenStr = "gnuplot -persistent";
	const char * commandOpen = commandOpenStr.c_str();
	FILE * gnuplotPipe = popen (commandOpen,"w");
	string commandStrPic1 = "set term png size 1600,800";
    string commandStrPic2 = "set output './pics/sphaleron.png'";
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
	
struct paramsStruct {double Y_0;};
struct paramsStruct p;
	
//Energy integrand
double EIntegrand (double x, void * parameters)
	{
	if (x<0.0)
		{
		printf ("error in EIntegrand, R<0");
		return 0.0;
		}
	else
		{
		struct paramsStruct * params = (struct paramsStruct *)parameters;
		double Y_0 = (params->Y_0);
		double y_R[4] = { Y_0, 0.0, 1.0, 0.0};
		int status;
		double t = 1.0e-16;
		gsl_odeiv2_system syse = {func, jac, 4, &paramsVoid};
		gsl_odeiv2_driver * de = gsl_odeiv2_driver_alloc_yp_new (&syse,  gsl_odeiv2_step_rk8pd, 1.0e-9, 1.0e-9, 0.0);
		status = gsl_odeiv2_driver_apply (de, &t, x, y_R);
		if (status != GSL_SUCCESS)
			{
			printf ("error in EIntegrand, return value=%d\n", status);
			gsl_odeiv2_driver_free (de);
			return (double)status;
			}
		else
			{
			double to_return = 4.0*pi*pow(x,2.0)*( 0.5*pow(y_R[1],2.0) + 0.5*pow(y_R[0],2.0) - 0.25*pow(y_R[0],4.0) );
			gsl_odeiv2_driver_free (de);
			return to_return;
			}
		}
	}

int main()
{
/* first of all just writing something to solve the ode as an initial value plot
only once that works will i add the newton raphson loop to do the shooting
to do:
- define gsl functions for first order odes
- define jacobian functions
- define boundary values
- follow gsl example for implementation

next, have to implement shooting by applying newton's method to F, the difference between the final boundary value and the one shot at
- define function F
- write loop with convergence test
- write newton's method
- print tests for convergence on each iteration

then,
- calculate energy of solution
- calculate number of particles
- find negative mode
*/
paramsVoid = {};

gsl_odeiv2_system sys = {func, jac, 4, &paramsVoid};

double F = 1.0, dF;
double aim = 0.0;
double closeness = 1.0e-9;
double t0 = 1.0e-16, t1 = 10.0;
const unsigned int steps = 1e3;
double h = t1-t0;
h /= (double)steps;
vector<double> y0Vec(steps+1), y2Vec(steps+1);
unsigned int runsCount = 0;

double Y0 = 4.337;//initial guess
//cout << "initial y0: ";
//cin >> Y0;

printf("%16s%16s%16s%16s%16s%16s%16s%16s\n","run","steps","y(t1)","yMin","F-aim","Y0Old","Y0New","-F/dF");

while (absolute(F-aim)>closeness)
	{
	runsCount++;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_yp_new (&sys,  gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

	double t = t0, ti;
	double y0[4] = { Y0, 0.0, 1.0, 0.0};
	double y[4];
	memcpy(y, y0, sizeof(y0));
	double yMin = absolute(y[0]-aim);
	y0Vec[0] = y[0];
	y2Vec[0] = y[2];
	int status;
	unsigned int i, iMin;
	
	for (i = 1; i <= steps; i++)
		{
		ti =(double)i;
		ti *= h;
		ti += t0;
		status = gsl_odeiv2_driver_apply (d, &t, ti, y);
		if (status != GSL_SUCCESS)
			{
			printf ("error, return value=%d\n", status);
			printf ("i = %3i, t = %3g\n",i,t);
			break;
			}
		if ((y[0]-aim)<yMin && (y[0]-aim)>0.0)
			{
			yMin = y[0]-aim;
			iMin = i;
			}
		y0Vec[i] = y[0];
		y2Vec[i] = y[2];
		//printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
		if ((y[0]-aim)<0.0)
			{
			iMin = i;
			if ((y[0]-aim)<-0.2) { break; }
			}
		}
	if (status != GSL_SUCCESS){break;}
		
	F = y0Vec[iMin]-aim; //as final boundary condition is y(t1)=0.0;
	dF = y2Vec[iMin];
	printf("%16i%16i%16g%16g%16g%16.12g",runsCount,i,y[0],yMin,F-aim,Y0);
	if (absolute(dF)>2.0e-16)
		{
		Y0 += -F/dF;
		printf("%16.12g%16g\n",Y0,-F/dF);
		}
	gsl_odeiv2_driver_free (d);
	if (i==(steps+1))
		{
		F = y[0];
		}
	}
	
p  = {Y0};

//finding E
double E, Eerror;
gsl_function E_integrand;
E_integrand.function = &EIntegrand;
E_integrand.params = &p;
gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
gsl_integration_qag(&E_integrand, 1.0e-16, 10.0, 1.0e-10, 1.0e-9, 1e4, 4, w, &E, &Eerror);
gsl_integration_workspace_free(w);
if (Eerror>1.0e-8) { cout << "E error = " << Eerror << endl;}
else { cout << "E = " << E << endl;}
	
string filename = "./data/sphaleron.dat";
simplePrintVector(filename,y0Vec);
printf("\n%16s%20s%20s\n","Solution printed: ",filename.c_str()," ./pics/sphaleron.png");
gpSimple(filename);

return 0;
}
