// error classes and function declarations

#ifndef __ERROR_H_INCLUDED__
#define __ERROR_H_INCLUDED__

#include <string>
#include <iostream>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
error handling
-------------------------------------------------------------------------------------------------------------------------*/

class SimpleError {
public:
	virtual string message() const = 0;
	virtual ~SimpleError();
};

ostream& operator<<(ostream& s, const SimpleError& e);

#endif // __ERROR_H_INCLUDED__
