// declarations of some very simple functions and classes

#ifndef __SIMPLE_H_INCLUDED__
#define __SIMPLE_H_INCLUDED__

#include <string>

using namespace std;

template <class T>
string numberToString ( T Number );

template <class T>
T stringToNumber ( const string& Text );

#endif // __SIMPLE_H_INCLUDED__
