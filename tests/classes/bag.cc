//program to test new classes

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "simple.h"
#include "bag.h"

/*-------------------------------------------------------------------------------------------------------------------------
fallible
-------------------------------------------------------------------------------------------------------------------------*/

template <class T>
Fallible<T>::UsedInInvalidState::UsedInInvalidState() {}

template <class T>
string Fallible<T>::UsedInInvalidState::message() const {
	return "Fallible object used in invalid state";
}

template <class T>
inline Fallible<T>::operator T() const {
if (failed()) throw UsedInInvalidState();
return instance;
}

/*-------------------------------------------------------------------------------------------------------------------------
bag errors
-------------------------------------------------------------------------------------------------------------------------*/

string BagError::NotFound::message() const{
	return "Parameter " + pName + " not found.";
}

/*-------------------------------------------------------------------------------------------------------------------------
parameter
-------------------------------------------------------------------------------------------------------------------------*/

template<class T>
Parameter<T>::Parameter() : name(), value() {}

template<class T>
Parameter<T>::Parameter(const string& pName, const T& pValue) : name(pName), value(pValue) {}

template<class T>
void Parameter<T>::copy(const Parameter& p) {
	name = p.name;
	value = p.value;
}

template<class T>
Parameter<T>::Parameter(const Parameter& p) {
	copy(p);
}

template<class T>
Parameter<T>& Parameter<T>::operator=(const Parameter& rhs) {
	copy(rhs);
	return *this;
}


template<class T>
bool Parameter<T>::empty() const {
return name.empty();
}

template<class T>
istream& operator>>(istream& is, Parameter<T>& p) {
	is >> p.name >> p.value;
	return is;
}

template<class T>
ostream& operator<<(ostream& os,const Parameter<T>& p) {
	return os << left << setw(20) << p.name << setw(20) << p.value << endl;
}

/*-------------------------------------------------------------------------------------------------------------------------
parameter bag
-------------------------------------------------------------------------------------------------------------------------*/

template <class T>
ParameterBag<T>::ParameterBag(): numParams(0), parameters() {}

template <class T>
void ParameterBag<T>::empty() {
	numParams = 0;
	parameters.clear();
}

template <class T>
void ParameterBag<T>::copy(const ParameterBag<T>& b) {
	empty();
	parameters = b.parameters;
	numParams = parameters.size();
}

template <class T>
ParameterBag<T>::ParameterBag(const ParameterBag<T>& b) {
	copy(b);
}

template <class T>
ParameterBag<T>& ParameterBag<T>::operator=(const ParameterBag<T>& b) {
	copy(b);
	return *this;
}

template<class T>
uint ParameterBag<T>::size() const {
	return numParams;
}

template<class T>
void ParameterBag<T>::add(const Parameter<T>& p) {
	parameters.push_back(p);
	numParams++;
}

template<class T>
Fallible <Parameter<T>* > ParameterBag<T>::find(const string& s) {
	uint it = 0;
	while(it<size()) {
		Parameter<T>& p = parameters[it];
		if(s.compare(p.name)==0) {
			return &p;
		}
		it++;
	}
	return Fallible< Parameter<T>* >();
}

template<class T>
Fallible <Parameter<T> > ParameterBag<T>::find(const string& s) const {
	uint it = 0;
	while(it<size()) {
		Parameter<T> p = parameters[it];
		if(s.compare(p.name)==0) {
			return p;
		}
		it++;
	}
	return Fallible< Parameter<T> >();
}

template <class T>
void ParameterBag<T>::set(const Parameter<T>& p) {
	Fallible <Parameter <T>* > f = find(p.name);
	if (f.valid()) {
		Parameter <T>* g = f;
		*g = p;
	}
	else add(p);
}

template <class T>
T ParameterBag<T>::operator()(const string& pName) const{
	Fallible <Parameter <T> > f = find(pName);
	if (f.valid()) {
		Parameter<T> p = f;
		return p.value;
	}
	BagError::NotFound e(pName);
	throw e;
}

template <class T>
void ParameterBag<T>::reset() {
	empty();
}

template<class T>
istream& operator>>(istream& is, ParameterBag<T>& b) {
	b.reset();
	uint it = 0;
	Parameter<T> p;
	while (!is.eof()) {
		is >> p;
		b.set(p);
		it++;
	}
	return is;
}

template<class T>
ostream& operator<<(ostream& os,const ParameterBag<T>& b) {
	if (b.size()>0) {
		Parameter<T> p;
		for(uint it=0; it<b.size(); it++) {
			if (!b.parameters[it].empty()) {
				p = b.parameters[it];
				os << p;
			}
		}
	}
	return os;
}

/*-------------------------------------------------------------------------------------------------------------------------
explicit template instantiation
-------------------------------------------------------------------------------------------------------------------------*/

template class Parameter<uint>;
template class Parameter<double>;
template ostream& operator<< <uint>(ostream&, const Parameter<uint>&);
template istream& operator>> <uint>(istream&, Parameter<uint>&);
template ostream& operator<< <double>(ostream&, const Parameter<double>&);
template istream& operator>> <double>(istream&, Parameter<double>&);

template class ParameterBag<uint>;
template class ParameterBag<double>;
template ostream& operator<< <uint>(ostream&, const ParameterBag<uint>&);
template istream& operator>> <uint>(istream&, ParameterBag<uint>&);
template ostream& operator<< <double>(ostream&, const ParameterBag<double>&);
template istream& operator>> <double>(istream&, ParameterBag<double>&);

