/*
	program to test thetaT.cc and thetaT.h
*/

#include <iostream>
#include <complex>
#include "potentials.h"
#include "thetaT.h"

using namespace std;

typedef complex<double> comp;

int main() {

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
