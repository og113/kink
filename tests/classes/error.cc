//defintions of error classes and functions

#include <string>
#include "error.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
error handling
-------------------------------------------------------------------------------------------------------------------------*/

extern SimpleError::~SimpleError(){}

extern ostream& operator<<(ostream&s, const SimpleError& e) {
	return s << e.message() << endl;
}

extern BagError::NotFound::NotFound(const string& s) {
	pName = s;
}

extern string BagError::NotFound::message() const{
	return "Parameter " + pName + " not found.";
}
