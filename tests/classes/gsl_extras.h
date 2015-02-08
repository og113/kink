/*
	declarations for functions using gnu standard library functions
*/

#ifndef __GSL_EXTRAS_H_INCLUDED__
#define __GSL_EXTRAS_H_INCLUDED__

#include <complex>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
generic functions packaging gsl functions
	- rootFinder
	- brentRootFinder
	- brentMinimum
	- rootsFn
-------------------------------------------------------------------------------------------------------------------------*/

//function to find a root of function FDF, given initial guess, using newton method
double rootFinder(gsl_function_fdf * xFDF, double rootGuess);
	
//function to find a root of function F, given initial guess, and lower and upper bounds, using brent method
double brentRootFinder(gsl_function * xF, const double & rootGuess, const double & rootLower, const double &rootUpper);

//function to find the minimum of a gsl function F, given initial guess, and lower and upper bounds, using brent method
double brentMinimum (gsl_function * xF, const double & minimumGuess, const double & minimumLower, const double & minimumUpper);

//function to give the three roots of FDF given lower and upper limits on them and a number of loops to try
vector <double> rootsFn (gsl_function_fdf * xFDF, const double & lowLimit, const double & highLimit,\
							const unsigned int & rootLoops);
	
#endif // __GSL_EXTRAS_H_INCLUDED__
