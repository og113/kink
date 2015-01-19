//defintions of error classes and functions

#include <string>
#include "error.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
error handling
-------------------------------------------------------------------------------------------------------------------------*/

SimpleError::~SimpleError(){}

ostream& operator<<(ostream&s, const SimpleError& e) {
	return s << e.message() << endl;
}

