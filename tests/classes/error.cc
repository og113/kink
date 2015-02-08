//defintions of error classes and functions

#include <string>
#include "error.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
base class for error handling
-------------------------------------------------------------------------------------------------------------------------*/

SimpleError::~SimpleError(){}

ostream& operator<<(ostream&s, const SimpleError& e) {
	return s << e.message() << endl;
}

/*-------------------------------------------------------------------------------------------------------------------------
file related errors
-------------------------------------------------------------------------------------------------------------------------*/

string FileError::StreamNotGood::message() const{
	return "Stream for " + filename + " not good.";
}
