// definitions of some very simple functions and classes

#include <string>
#include <sstream>
#include <cmath>
#include "simple.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
generic functions
-------------------------------------------------------------------------------------------------------------------------*/

//to convert number to string, usage is string str = NumberToString<number type>(x);
template <class T>
string numberToString ( const T& Number )
	{
	stringstream ss;
	ss << Number;
	return ss.str();
	}

//to convert string to number, usage is (number type) x = StringToNumber<number type>(str);
template <class T>
T stringToNumber ( const string& Text )
	{
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
	}
	
//gives absolute measure of difference between two numbers
double absDiff (const double& numA, const double& numB) {
	return 2.0*abs(numA-numB)/sqrt(numA*numA+numB*numB);
}

/*-------------------------------------------------------------------------------------------------------------------------
explicit template instantiation
-------------------------------------------------------------------------------------------------------------------------*/

template string numberToString<int>(const int&);
template string numberToString<unsigned int>(const unsigned int&);
template string numberToString<double>(const double&);

template int stringToNumber<int>(const string&);
template unsigned int stringToNumber<unsigned int>(const string&);
template double stringToNumber<double>(const string&);

