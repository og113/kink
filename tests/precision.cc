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

using namespace std;

int main()
{

double a = 1.12387623476234;

cout << a << endl;
cout.precision(16);
cout << a << endl;
cout << a << endl;
cout << setprecision(5) << a << endl;
cout << a << endl;

fstream f;
f.open("tests/precision.txt", ios::out);
f.precision(10);
f << a << endl;
f << a << endl;
f.close();
	
return 0;
}
