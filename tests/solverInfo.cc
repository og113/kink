//program to check how the info of Eigen solvers works
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;

int main()
{
vec x(4);
vec b(4);
b << 1, 2, 3, 4;

spMat m(4,4);
m.insert(0,0) = 1;
m.insert(1,1) = 1;
m.insert(2,2) = 1;
m.insert(3,3) = 2;
m.makeCompressed();

Eigen::SparseLU<spMat> solver;
solver.analyzePattern(m);
if(solver.info()!=Eigen::Success) {
cout << "pattern analysis failed" << endl;
return 0;
}

solver.factorize(m);
if(solver.info()!=Eigen::Success) {
cout << "factorization failed" << endl;
return 0;
}

x = solver.solve(b);
if(solver.info()!=Eigen::Success) {
cout << "solving failed" << endl;
return 0;
}

cout << solver.info() << endl;

cout << x << endl;

return 0;
}
