/*
	declarations for some extra functions dealing with std::map
*/

#ifndef __MAP_EXTRAS_H_INCLUDED__
#define __MAP_EXTRAS_H_INCLUDED__

#include <fstream>
#include <iostream>
#include <string>
#include <map>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
<< and >> for pair and map
-------------------------------------------------------------------------------------------------------------------------*/

template <class Value>
ostream& operator<<(ostream&, const map<string,Value>&);

template <class Value>
ifstream& operator>>(ifstream&, map<string,Value>&);

#endif // __MAP_EXTRAS_H_INCLUDED__
