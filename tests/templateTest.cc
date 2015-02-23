//quick test example function to see whether I can pass a function of type T
#include <stdio.h>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

using namespace std;

typedef complex<double> comp;

double quickFn (const double & x) {
	return x*x;
}

//params
double epsilon, A;
struct params_for_V {double epsi; double aa;};
struct params_for_V paramsV {epsilon, A};

//V1
template <class T> T V1 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	return pow(pow(phi,2)-1.0,2.0)/8.0 - epsi*(phi-1.0)/2.0;
	}	

//Z, for V2
template <class T> T Z (const T phi)
	{
	return exp(-pow(phi,2.0))*(phi + pow(phi,3.0) + pow(phi,5.0));
	}
	
//Y, for V2
template <class T> T Y (const T phi)
	{
	return 0.5*pow(phi+1.0,2.0);
	}
	
//Z, for V2
template <class T> T X (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	double aa = (params->aa);
	printf("%8g\n\n",aa);
	printf("%8g\n\n",epsi);
	return (phi-1.0)/aa;
	}
	
//V2
template <class T> T V2 (const T phi, void * parameters = &paramsV)
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	double epsi = (params->epsi);
	double aa = (params->aa);
	return 0.5*pow(phi+1.0,2.0)*(1.0-epsi*Z((phi-1.0)/aa));
	}
comp V2c (const comp phi) { return V2(phi); }

double (*fnp) (const double x, void * parameters);
comp (*fnpc) (const comp x, void * parameters);

int
main (void)
{
  epsilon = 0.75;
  A = 0.4;
  paramsV  = {epsilon, A};
  //printf("%8g\n\n",A);
  //printf("%8g\n\n",epsilon);
  //printf("%8g\n\n",Z(-1.0));
  //printf("%8g\n\n",Y(-1.0));
  //printf("%8g\n\n",X(-1.0));
  //printf("%8g\n\n",V2(-1.0));
  
  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = 1.2, m_expected = 1.6;
  double a = 0.5, b = 5.0;
  gsl_function F;

  F.function = &V2;
  F.params = &paramsV;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

  printf ("using %s method\n",
          gsl_min_fminimizer_name (s));

  printf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "min",
          "err", "err(est)");

  printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
          iter, a, b,
          m, m - m_expected, b - a);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status 
        = gsl_min_test_interval (a, b, 0.001, 0.0);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] "
              "%.7f %+.7f %.7f\n",
              iter, a, b,
              m, m - m_expected, b - a);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free (s);
  
  template <class T>
  typedef T (Potential*) (const T&);
  Potential<double> pot = &quickFn;
  cout << pot(9) << endl;
  

  return status;
}
