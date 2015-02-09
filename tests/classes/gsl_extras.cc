/*
	definitions for functions using gnu standard library functions
	must compile with -lgsl and -lgslcblas
*/

#include <complex>
#include <iostream>
#include <vector>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "gsl_extras.h"
#include "potentials.h"

/*-------------------------------------------------------------------------------------------------------------------------
generic functions packaging gsl functions
	- double rootFinder(gsl_function_fdf * xFDF, double rootGuess)
	- double brentRootFinder(gsl_function * xF, const double & rootGuess, const double & rootLower, const double &rootUpper)
	- double brentMinimum (gsl_function * xF, const double & minimumGuess, const double & minimumLower, const double & minimumUpper)
	- vector <double> rootsFn (gsl_function_fdf * xFDF, const double & lowLimit, const double & highLimit,\
							const unsigned int & rootLoops)
-------------------------------------------------------------------------------------------------------------------------*/

//function to find a root of function FDF, given initial guess, using newton method
double rootFinder(gsl_function_fdf * xFDF, double rootGuess)
	{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	double x0, x = rootGuess;

	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc (T);
	gsl_root_fdfsolver_set (s, xFDF, x);

	do
		{
		  iter++;
		  status = gsl_root_fdfsolver_iterate (s);
		  x0 = x;
		  x = gsl_root_fdfsolver_root (s);
		  status = gsl_root_test_delta (x, x0, 0, DBL_MIN);
		}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fdfsolver_free (s);
	return x;
	}
	
//function to find a root of function FDF, given initial guess, and lower and upper bounds, using brent method
double brentRootFinder(gsl_function * xF, const double & rootGuess, const double & rootLower, const double &rootUpper)
	{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double x = rootGuess;
	double x_lo = rootLower;
	double x_hi = rootUpper;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, xF, x_lo, x_hi);

	do
		{
		  iter++;
		  status = gsl_root_fsolver_iterate (s);
		  x = gsl_root_fsolver_root (s);
		  x_lo = gsl_root_fsolver_x_lower (s);
		  x_hi = gsl_root_fsolver_x_upper (s);
      	  status = gsl_root_test_interval (x_lo, x_hi,
                                       0, DBL_MIN);
		}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
	return x;
	}
	
//function to find the minimum of a gsl function F
double brentMinimum (gsl_function * xF, const double & minimumGuess, const double & minimumLower, const double & minimumUpper)
	{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m = minimumGuess;
	double tempLower, tempUpper;

	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc (T);
	gsl_min_fminimizer_set (s, xF, m, minimumLower, minimumUpper);
	
	do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      tempLower = gsl_min_fminimizer_x_lower (s);
      tempUpper = gsl_min_fminimizer_x_upper (s);

      status = gsl_min_test_interval (tempLower, tempUpper, 1.0e-12, DBL_MIN);

      if (status == GSL_SUCCESS)
      	{
        //printf ("brentMinimum converged\n");
        }
    }
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_min_fminimizer_free (s);
	
	return m;
	}

//function to give the three roots of FDF given lower and upper limits on them and a number of loops to try
vector <double> rootsFn (gsl_function_fdf * xFDF, const double & lowLimit, const double & highLimit,\
							const unsigned int & rootLoops)
	{
	vector <double> roots;
	for (unsigned int j=0;j<rootLoops;j++)
		{
		double x = lowLimit+(abs(highLimit)+abs(lowLimit))*j/(rootLoops-1.0);
		x = rootFinder(xFDF,x);
	
		if (j==0) roots.push_back(x);
		else
			{
			unsigned int test = 0;
			for (unsigned int k=0;k<roots.size();k++) if (abs(x-roots[k])>1.0e-12) test++;
			if (test==roots.size()) roots.push_back(x);
			}
		}
	
		if (roots.size()!=3) cout << "rootsFn error: only found " << roots.size() << " roots, not 3" << endl;
	return roots;
	}

