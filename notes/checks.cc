/*
	definitions for checks class
*/

#include <iostream>
#include "checks.h"

bool Check::operator()(const double& test) {
	if (test<closeness) return true;
	else {
		cerr << message << endl;
		cerr << "test = " test << " > " << closeness << endl;
		return false;
	}
}
