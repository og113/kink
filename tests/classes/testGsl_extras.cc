/*
	main for program to test gsl_extras.cc and gsl_extras.h
*/

#include <iostream>
#include "potentials.h"
#include "gsl_extras.h"

using namespace std;

int main() {
cout << "test: " << endl;

struct params_for_V paramsV;
paramsV.epsi = 0.05;
paramsV.aa = 0.4;

struct void_params paramsVoid;

double d = 0.9;
comp c(-1.0,0.0);

PotentialType pot = &V2<comp>;

Potential classV(paramsV,pot);

cout << classV(c) << endl;

return 0;
}
