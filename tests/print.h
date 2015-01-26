/*
declarations for printing functions

conventions:
	- string filename should always be first parameter
*/

#ifndef __PRINT_H_INCLUDED__
#define __PRINT_H_INCLUDED__

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "error.h"

typedef unsigned int uint;
typedef Eigen::VectorXd vec;

/*-------------------------------------------------------------------------------------------------------------------------
struct holding print options
-------------------------------------------------------------------------------------------------------------------------*/
struct PrintOptions {
	
};
/*-------------------------------------------------------------------------------------------------------------------------
error handling for print functions
-------------------------------------------------------------------------------------------------------------------------*/

class PrintError {
public:
	class FileUnopen: public SimpleError{
	public:
		FileUnopen(const string& s) : filename(s) {}		// constructor
		virtual string		message() const;				// message to be passed for printing
	private:
		string	filename;									// name of file unopened
	};
	class EmptyVector: public SimpleError{
	public:
		FileUnopen() {}										// constructor
		virtual string		message() const;				// message to be passed for printing
	}
};

/*-------------------------------------------------------------------------------------------------------------------------
simple print functions for vectors
-------------------------------------------------------------------------------------------------------------------------*/

typedef <class Vector>
void print (const string& filename, const Vector& vector, const struct PrintOptions& options);

#endif // __PRINT_H_INCLUDED__
