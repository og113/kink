/*
	a quick program to perform the integrations needed to get E(T) for the low energy periodic instantons of the upside-down phi^4 potential
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

using namespace std;

#define pi    3.1415926535898
#define euler 0.5772156649015

struct paramStruct {double t;};

double omega(double k) {
	return sqrt(1.0+pow(k,2.0));
}

double F_integrand(double k, void* pS) {
	struct paramStruct* params = (struct paramStruct *)pS;
	double t = (params->t);
	return 8.0*pow(k,2.0)/omega(k)/(exp(omega(k)*t)-1.0);
}

double F(const double& t) {
	struct paramStruct pS;
	pS.t = t;
	double F_error, F_result;
	gsl_function F_integrand_gsl;
	F_integrand_gsl.function = &F_integrand;
	F_integrand_gsl.params = &pS;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
	gsl_integration_qag(&F_integrand_gsl, 0.0, 1.0e2, 1.0e-16, 1.0e-8, 1e4, 4, w, &F_result, &F_error);
	gsl_integration_workspace_free(w);
	if (F_error>1.0e-8) { cout << "F error = " << F_error << endl;}
	return F_result;
}

double rho_sqrd(const double& t) {
	return 4.0*exp(-2.0*euler-2.0-F(t));
}

double E_integrand(double k, void* pS) {
	struct paramStruct* params = (struct paramStruct *)pS;
	double t = (params->t);
	return pow(k,2.0)*exp(omega(k)*t)/pow((exp(omega(k)*t)-1.0),2.0);
}

typedef double (functionType*)(double,void*);

double integrate(double t, ) {
	struct paramStruct pS;
	pS.t = t;
	double E_error, E_result;
	gsl_function E_integrand_gsl;
	E_integrand_gsl.function = &E_integrand;
	E_integrand_gsl.params = &pS;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
	gsl_integration_qag(&E_integrand_gsl, 0.0, 1.0e2, 1.0e-16, 1.0e-8, 1e4, 4, w, &E_result, &E_error);
	gsl_integration_workspace_free(w);
	if (E_error>1.0e-8) { cout << "E error = " << E_error << endl;}
	return 32.0*pow(pi,2.0)*rho_sqrd(t)*E_result;
}

double S(double t) {
	return (8.0*pow(pi,2.0)/3.0)*(1.0-(1.0/6.0)*rho_sqrd(t)*(2*log(pow(rho_sqrd(t),0.5)/2.0)+2.0*euler+1+F(t)));
}

int main() {

double T = 1.5;

cout << left;
cout << setw(20) << "T" << setw(20) << "F" << setw(20) << "rho_sqrd" << setw(20) << "E" << setw(20) << "S" << setw(20) << "W" << endl;
for (unsigned int j=0; j<21; j++) {
	cout << setw(20) << T;
	cout << setw(20) << F(T);
	cout << setw(20) << rho_sqrd(T);
	cout << setw(20) << E(T)/18.9;
	cout << setw(20) << S(T)/(8.0*pow(pi,2.0)/3.0);
	cout << setw(20) << (S(T)-E(T)*T)/(8.0*pow(pi,2.0)/3.0) << endl;
	T += 0.1;
}

return 0;
}
