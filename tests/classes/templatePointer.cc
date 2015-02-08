/*
program to test interaction of templates and function pointers
*/

#include <iostream>
#include <string>
#include "templatePointer.h"

using namespace std;

template <class T>
T function(T x) {
	return x;
}

template double function<double>(double);
template int function<int>(int);
template string function<string>(string);
