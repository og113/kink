//parameters and functions for pi.cc

#ifndef __BAG_H_INCLUDED__
#define __BAG_H_INCLUDED__

#include <iostream>
#include <vector>
#include "error.h"

using namespace std;

typedef unsigned int uint;

/*-------------------------------------------------------------------------------------------------------------------------
fallible
	- useful return type for search functions
-------------------------------------------------------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------------------------------------------------------
bag errors
-------------------------------------------------------------------------------------------------------------------------*/

class BagError {
public:
	class NotFound: public SimpleError{
	public:
		NotFound(const string& s) : pName(s) {}		// constructor
		virtual string		message() const;		// message to be passed for printing
	private:
		string	pName;								// name of parameter not found
	};
};

/*-------------------------------------------------------------------------------------------------------------------------
parameter
-------------------------------------------------------------------------------------------------------------------------*/

//defines a parameter
template<class T>
class Parameter {
public:
	Parameter();										// empty constructor
	Parameter(const string& pName, const T& pValue);	// constructor
	Parameter(const Parameter&);						// copy constructor
	Parameter& operator=(const Parameter&);				// assign parameter
	~Parameter() {}										// destructor
	bool		empty() const;							// gives true if empty
	
	string name;										// name of parameter
	T value;											// value of parameter
	
private:
	void copy(const Parameter& p);						// private copy function
};

template<class T>
istream& operator>>(istream& is, Parameter<T>& p);

template<class T>
ostream& operator<<(ostream& os,const Parameter<T>& p);

/*-------------------------------------------------------------------------------------------------------------------------
parameter bag
-------------------------------------------------------------------------------------------------------------------------*/
template <class T>
class ParameterBag;

template <class T>
istream& operator>>(istream&, ParameterBag<T>&);
template <class T>
ostream& operator<<(ostream&, const ParameterBag<T>&);

template <class T>
class ParameterBag {
public:
	ParameterBag();											// empty constructor
	ParameterBag(const ParameterBag&);						// copy constructor
	ParameterBag& operator=(const ParameterBag&);			// assign parameter bag
	~ParameterBag() {}										// destructor
	T operator()(const string& pName) const;				// to get a named T parameter
	void				set(const Parameter<T>& p);			// set parameter
	void				reset();							// set to empty
	uint				size() const;						// gets numParams
	
	friend istream& operator>> <T>(istream&, ParameterBag&);
	friend ostream& operator<< <T>(ostream&, const ParameterBag&);
	
private:
	uint					numParams;						// number of parameters
	vector< Parameter<T> >	parameters;						// array of parameters
	void					copy(const ParameterBag& p);	// private copy function
	void					empty();						// empties bag
	void					add(const Parameter<T>& p);		// adds parameter
	Fallible< Parameter<T>* >find(const string& s);			// if present, returns parameter
	Fallible< Parameter<T> > find(const string& s) const;	// if present, returns parameter
};

template <class T>
istream& operator>>(istream& is, ParameterBag<T>& b);

template <class T>
ostream& operator<<(ostream& os, const ParameterBag<T>& b);

/*-------------------------------------------------------------------------------------------------------------------------
parameter bag separating primary and secondary parameters
-------------------------------------------------------------------------------------------------------------------------*/

class ParameterBag_ps: public ParameterBag<double>, ParameterBag<uint> {
public:

private:
};

#endif // __BAG_H_INCLUDED__
