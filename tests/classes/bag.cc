//program to test new classes

#ifndef __SIMPLE_H_INCLUDED__
#define __SIMPLE_H_INCLUDED__

#ifndef __BAG_H_INCLUDED__
#define __BAG_H_INCLUDED__

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
extern Fallible<T>::UsedInInvalidState::UsedInInvalidState() {}

template <class T>
extern string Fallible<T>::UsedInInvalidState::message() const {
	return "Fallible object used in invalid state";
}

template <class T>
inline
extern Fallible<T>::operator T() const {
if (failed()) throw UsedInInvalidState();
return instance;
}

/*-------------------------------------------------------------------------------------------------------------------------
parameter
-------------------------------------------------------------------------------------------------------------------------*/

template<class T>
extern Parameter<T>::Parameter() : name(), value() {}

template<class T>
extern Parameter<T>::Parameter(const string& pName, const T& pValue) : name(pName), value(pValue) {}

template<class T>
extern void Parameter<T>::copy(const Parameter& p) {
	name = p.name;
	value = p.value;
}

template<class T>
extern Parameter<T>::Parameter(const Parameter& p) {
	copy(p);
}

template<class T>
extern Parameter<T>& Parameter<T>::operator=(const Parameter& rhs) {
	copy(rhs);
	return *this;
}


template<class T>
extern bool Parameter<T>::empty() const {
return name.empty();
}

template<class T>
extern istream& operator>>(istream& is, Parameter<T>& p) {
	is >> p.name >> p.value;
	return is;
}

template<class T>
extern ostream& operator<<(ostream& os,const Parameter<T>& p) {
	return os << left << setw(20) << p.name << setw(20) << p.value << endl;
}

/*-------------------------------------------------------------------------------------------------------------------------
parameter bag
-------------------------------------------------------------------------------------------------------------------------*/

template <class T>
extern ParameterBag<T>::ParameterBag(): numParams(0), parameters() {}

template <class T>
extern void ParameterBag<T>::empty() {
	numParams = 0;
	parameters.clear();
}

template <class T>
extern void ParameterBag<T>::copy(const ParameterBag<T>& b) {
	empty();
	parameters = b.parameters;
	numParams = parameters.size();
}

template <class T>
extern ParameterBag<T>::ParameterBag(const ParameterBag<T>& b) {
	copy(b);
}

template <class T>
extern ParameterBag<T>& ParameterBag<T>::operator=(const ParameterBag<T>& b) {
	copy(b);
	return *this;
}

template<class T>
extern uint ParameterBag<T>::size() const {
	return numParams;
}

template<class T>
extern void ParameterBag<T>::add(const Parameter<T>& p) {
	parameters.push_back(p);
	numParams++;
}

template<class T>
extern Fallible <Parameter<T>* > ParameterBag<T>::find(const string& s) {
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
extern Fallible <Parameter<T> > ParameterBag<T>::find(const string& s) const {
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
extern void ParameterBag<T>::set(const Parameter<T>& p) {
	Fallible <Parameter <T>* > f = find(p.name);
	if (f.valid()) {
		Parameter <T>* g = f;
		*g = p;
	}
	else add(p);
}

template <class T>
extern ParameterBag::ParameterBag(const int argc&, const char** argv) {
	if (argc>1) {
		for (unsigned int j=0; j<(int)(argc/2); j++) {
			string pName = argv[(int)(argc/2) + 1];
			pName = pName.substr(1,pName.size());
			Parameter p(pName,stringToNumber(argv[(int)(argc/2)+2]));
			set(p);
		}
	}
}

template <class T>
extern T ParameterBag<T>::operator()(const string& pName) const{
	Fallible <Parameter <T> > f = find(pName);
	if (f.valid()) {
		Parameter<T> p = f;
		return p.value;
	}
	BagError::NotFound e(pName);
	throw e;
}

template <class T>
extern void ParameterBag<T>::reset() {
	empty();
}

template<class T>
extern istream& operator>>(istream& is, ParameterBag<T>& b) {
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
extern ostream& operator<<(ostream& os,const ParameterBag<T>& b) {
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
main
-------------------------------------------------------------------------------------------------------------------------*/

int main(){
Parameter<double> p("csi",1.0);
Parameter<double> q(p);
Parameter<double> r;
r = q;
cout << p << q << r;
Parameter<double> s("chi",2.0), t("pi",3.0), u("csi",0.0);
ParameterBag<double> b, c;
b.set(s);
b.set(t);
cout << b.size() << endl;
cout << b;
b.set(p);
b.set(u);
cout << b.size() << endl;
cout << b;
try{
cout << b("csi") << endl;
} catch (BagError::NotFound e) {
cerr << e;
}
try{cout << b("fi") << endl;
} catch (BagError::NotFound e) {
cerr << e;
}
b.reset();
b.set(r);
c = b;
cout << c;
c.reset();

ifstream file;
file.open("inputs");
file >> c;
file.close();
cout << c;

ofstream out;
out.open("outputs");
out << c;
out.close();


return 0;
}
