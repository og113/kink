// declarations of some very simple functions and classes

#ifndef __SIMPLE_H_INCLUDED__
#define __SIMPLE_H_INCLUDED__

#include <string>

using namespace std;

template <class T>
string numberToString ( const T& Number );

template <class T>
T stringToNumber ( const string& Text );

double absDiff (const double& numA, const double& numB);

#endif // __SIMPLE_H_INCLUDED__
