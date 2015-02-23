//program to initialize a vector and matrix to zero
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXcd vec;

int main()
{
vec v(3);
v = Eigen::VectorXcd::Zero(3);

mat m(3,3);
m  = Eigen::MatrixXd::Zero(3,3);

cout << v(0) << endl;

cout << m(0,0) << endl;

return 0;
}
