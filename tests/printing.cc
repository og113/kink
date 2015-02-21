//program to check fprintf
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
typedef Eigen::VectorXcd vec;

int main()
{
//primary parameters
unsigned int N = 80; //number of points in each spatial dimension
unsigned int Na = (int)(1.2*N);
unsigned int Nb = (int)(1.0*N);
unsigned int Nc = 2;
double R = 10.0; //size of bubble
double mass = 3.0; 
double lambda = 0.1;
double Tb = 1.2*R/2;
double angle = asin(Tb/R); //not a primary parameter, just used to make L
double L = 2*(1.5*Tb*tan(angle));
double theta = 0.0;

vec v(3);
v << 1.0, 2.0, 3.0;

printf("%12g",1.0);
printf("%12g\n",3.0);

	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","N","Na","Nb","Nc","L","Tb","R","mass","lambda","theta");
	printf("%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g\n",N,Na,Nb,Nc,L,Tb,R,mass,lambda,theta);

return 0;
}
