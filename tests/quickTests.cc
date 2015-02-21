/*
	quick program to do quick tests
*/

#include <iostream>
#include "files.h"

using namespace std;

template <class T>
T quickFn(const T& x) {
	return x*x;
}

double (*fnPtr)(const double& x);

typedef double (*DoublePotential)(const double&);

int main() {
fnPtr = &quickFn<double>;
DoublePotential pot = fnPtr;
cout << pot(26.0) << endl;

changeInputs("data/quickTestInputs","theta","0.01","inputs");

return 0;
}
