/*
	definitions for some extra functions dealing with std::map
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <map>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
<< and >> for map
-------------------------------------------------------------------------------------------------------------------------*/

template <class Value>
ostream& operator<<(ostream& os, const map<string,Value>& mp) {
	for (typename map<string,Value>::const_iterator it=mp.begin(); it!=mp.end(); ++it)
		os << left << setw(20) << (*it).first << setw(20) << (*it).second << endl;
	os << endl;
}

template <class Value>
ifstream& operator>>(ifstream& is, map<string,Value>& mp) {
	string line, temp;
	while(getline(is,line)) {
		if(line[0] == '#') continue;
		if(line.empty()) break;
		istringstream ss(line);
		ss >> temp;
		ss >> mp[temp];
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
explicit template instantiation
-------------------------------------------------------------------------------------------------------------------------*/

template ostream& operator<< <double>(ostream&, const map<string,double>&);
template ostream& operator<< <int>(ostream&, const map<string,int>&);
template ostream& operator<< <unsigned int>(ostream&, const map<string,unsigned int>&);
template ostream& operator<< <string>(ostream&, const map<string,string>&);

template ifstream& operator>> <double>(ifstream&, map<string,double>&);
template ifstream& operator>> <int>(ifstream&, map<string,int>&);
template ifstream& operator>> <unsigned int>(ifstream&, map<string,unsigned int>&);
template ifstream& operator>> <string>(ifstream&, map<string,string>&);
