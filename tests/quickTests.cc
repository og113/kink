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

unsigned int Na = 320;
unsigned int linearInt = (unsigned int)(Na/10);
cout << linearInt << endl;
double E_exact = 0.0;
for (unsigned int j=1; j<(linearInt+1); j++) E_exact += 0.1;
cout << E_exact << endl;
E_exact /= (double)linearInt;
cout << E_exact << endl;

return 0;
}
