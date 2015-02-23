//program to resize vector
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
typedef Eigen::VectorXd vec;

int main()
{
vec v(3);
v << 1, 2, 3;

v.conservativeResize(4);
v(3) = 0;

cout << v << endl;

return 0;
}
