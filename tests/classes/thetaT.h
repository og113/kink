/*
	declarations for functions specifically useful to the theta/T method
*/

#ifndef __THETAT_H_INCLUDED__
#define __THETAT_H_INCLUDED__

#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
specific functions F, FDF etc
	- VdV
	- dVddV
	- struct ec_params
	- ec (energy change: V(minima[1])-V(minima[0])-dE)
	- S1 integrand
	- rho integrand
-------------------------------------------------------------------------------------------------------------------------*/
//V FDF gsl function	
void VdV (double x, void * parameters, double * f, double* df);

//dV FDF gsl function	
void dVddV (double x, void * parameters, double * f, double* df);
	
//energy change parameter struct
struct ec_params {double aa; double minima0; double minima1; double de; };

//energy change F gsl function : V(minima[1])-V(minima[0])-dE
double ec (double epsi, void * parameters);
	
//S1 integrand
double s1Integrand (double x, void * parameters);

//rho integrand
double rhoIntegrand (double x, void * parameters);

/*-------------------------------------------------------------------------------------------------------------------------
	functions to find secondary parameters
		- epsilonFn
-------------------------------------------------------------------------------------------------------------------------*/

//program to find epsilon given gsls function df and dE
void epsilonFn (gsl_function * xF, gsl_function * xEC, double * xdE, double * xEpsilon, vector<double>* xMinima);
	
#endif // __THETAT_H_INCLUDED__
