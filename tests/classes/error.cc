//defintions of error classes and functions

#include <string>
#include "simple.h"
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
	return "Stream for " + Filename + " not good.";
}

/*-------------------------------------------------------------------------------------------------------------------------
index related errors
-------------------------------------------------------------------------------------------------------------------------*/

string IndexError::OutOfBounds::message() const{
	return "Index error: Out of Bounds. " + numberToString<unsigned int>(Index) + " not in (" \
		 + numberToString<unsigned int>(Min) + "," + numberToString<unsigned int>(Max) + ").";
}
