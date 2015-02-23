//program to check fprintf
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
int a = 1;
double b = 2.0;
complex<double> c(3.0,4.0);
complex<double> d;
double e;
complex<int> f(1,0);

e = a;

d = b*c + a*b;
e = a + b*a;

cout << d << endl;
cout << e << endl;

return 0;
}
