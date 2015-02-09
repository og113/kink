// error classes and function declarations

#ifndef __ERROR_H_INCLUDED__
#define __ERROR_H_INCLUDED__

#include <string>
#include <iostream>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
base class for error handling
-------------------------------------------------------------------------------------------------------------------------*/

class SimpleError {
public:
	virtual string message() const = 0;
	virtual ~SimpleError();
};

ostream& operator<<(ostream& s, const SimpleError& e);

/*-------------------------------------------------------------------------------------------------------------------------
file related errors
-------------------------------------------------------------------------------------------------------------------------*/

class FileError {
public:
	class StreamNotGood: public SimpleError{
	public:
		StreamNotGood(const string& s) : filename(s) {}		// constructor
		virtual string		message() const;			// message to be passed for printing
	private:
		string	filename;
	};
};

#endif // __ERROR_H_INCLUDED__
