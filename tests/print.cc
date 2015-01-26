/*
defintions for printing functions
*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "print.h"

/*-------------------------------------------------------------------------------------------------------------------------
error handling for print functions
-------------------------------------------------------------------------------------------------------------------------*/

string PrintError::FileUnopen::message() const{
	return "Print error: " + filename + " not opened properly.";
}

string PrintError::EmptyVector::message() const{
	return "Print error: empty vector passed.";
}

/*-------------------------------------------------------------------------------------------------------------------------
simple print functions for vectors
-------------------------------------------------------------------------------------------------------------------------*/

typedef <class T>
void print (const string& filename, const vector<T>& vector, const struct PrintOptions& options) {
	unsigned int length = vector.size();
	if (length==0) {
		PrintError::EmptyVector err;
		throw err;
	}
	ofstream os;
	os.open((filename).c_str());
	if (!os.good()) {
		PrintError::FileUnopen err(filename);
		throw err;
	}
	os.precision(16);
	os << left;
	for (unsigned int j=0; j<length; j++) os << setw(25) << vector[j] << endl;
	os.close();
}

/*-------------------------------------------------------------------------------------------------------------------------
explicit template instantiation
-------------------------------------------------------------------------------------------------------------------------*/

template void print (const string& filename, const vector<double>& vector, const struct PrintOptions& options);
template void print (const string& filename, const vector<int>& vector, const struct PrintOptions& options);
template void print (const string& filename, const vector<uint>& vector, const struct PrintOptions& options);
template void print (const string& filename, const vec& vector, const struct PrintOptions& options);

