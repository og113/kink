/*
	declarations for checks class
*/

#ifndef __CHECKS_H_INCLUDED__
#define __CHECKS_H_INCLUDED__

#include <string>
#include <iostream>

using namespace std;

class Check {
public:
	Check(const double & c, const string& s): closeness(c), message(s) {}
	~Check();
	bool operator()(const double& test);
private:
	double closeness;
	string message;
};


#endif // __CHECKS_H_INCLUDED__
