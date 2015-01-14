// error classes and function declarations
#include <string>
#include <iostream>

using namespace std;

typedef unsigned int uint;

/*-------------------------------------------------------------------------------------------------------------------------
error handling
-------------------------------------------------------------------------------------------------------------------------*/

class SimpleError {
public:
	virtual string message() const = 0;
	virtual ~SimpleError();
};

ostream& operator<<(ostream& s, const SimpleError& e);

//fallible template class as return type for search functions
template <class T>
class Fallible {
public:
	Fallible()			: isValid(false), instance() {}			// empty constructor - initialized to false
	Fallible(const T& t): isValid(true), instance(t) {}			// constructor - initialized to true
	bool failed() const { return !isValid;			  }			// true if invalid
	bool valid()  const { return isValid;			  }			// true if valid
	void invalidate()   { isValid = false;			  }			// make invalid
	operator T() const;											// provides conversion to T
	T elseDefaultTo(const T& defaultValue) const;				// value if valid, else defaultValue
	
	class UsedInInvalidState: public SimpleError {
	public:
		UsedInInvalidState();
		virtual string message() const;
	};
	
private:
	bool	isValid;
	T 		instance;
};

