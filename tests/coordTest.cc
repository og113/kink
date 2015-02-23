//program to check intCoord and coord
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include "pf.h"

using namespace std;

int main()
{

for (unsigned int j=0; j<3*4; j++)
	{
	unsigned int t = intCoord(j,0,NT);
	unsigned int x = intCoord(j,1,NT);
	cout << t << "    " << x << endl;
	}
cout << endl << endl;

for (unsigned int j=0; j<3*4; j++)
	{
	unsigned int t = intCoord(j,0,NT);
	double temp = (double)t;
	temp -= (double)Na;
	complex<double> t0 = b*temp + i*Tb;
	complex<double> x = coordA(j,1);
	cout << t << "     " << temp << "     " << t0 << "    " << x << endl;
	}
	
cout << "b    " << b << endl;
cout << "Tb    " << Tb << endl;
cout << "NT    " << NT << endl;
cout << "Na    " << Na << endl;
	
return 0;
}
