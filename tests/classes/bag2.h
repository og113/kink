//parameters and functions for pi.cc

#ifndef __BAG_H_INCLUDED__
#define __BAG_H_INCLUDED__

#include <iostream>
#include <vector>
#include <map>
#include "error.h"

using namespace std;

typedef unsigned int uint;

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
	- pair<string,T>
-------------------------------------------------------------------------------------------------------------------------*/

template<class Key, class Value>
istream& operator>>(istream& is, pair<Key,Value>& p);

template<class T>
ostream& operator<<(ostream& os,const pair<Key,Value>& p);

/*-------------------------------------------------------------------------------------------------------------------------
parameter bag
-------------------------------------------------------------------------------------------------------------------------*/
class ParameterBag;

istream& operator>>(istream&, ParameterBag&);
ostream& operator<<(ostream&, const ParameterBag&);

class ParameterBag {
public:
	ParameterBag();											// empty constructor
	ParameterBag(const ParameterBag&);						// copy constructor
	ParameterBag& operator=(const ParameterBag&);			// assign parameter bag
	~ParameterBag() {}										// destructor
	T operator[](const string& pName) const;				// to get a named T parameter
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
