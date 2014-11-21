//parameters and functions for pi.cc
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <ctype.h>
#include <limits> //for DBL_MIN etc
#include <cstring> //for memcpy
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "gnuplot_i.hpp"

using namespace std;

typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;
typedef Eigen::VectorXd vec;

struct void_params {};
struct void_params paramsVoid;

#define pi 3.14159265359

int func (double t, const double y[], double f[], void *params)
	{
	f[0] = y[1];
	f[1] = -2.0*y[1]/t + y[0] - pow(y[0],3.0);
	f[2] = y[3];
	f[3] = -2.0*y[3]/t + y[2] - 3.0*y[2]*pow(y[0],2.0);
	return GSL_SUCCESS;
	}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
	{
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 0, 2, 0.0);
	gsl_matrix_set (m, 0, 3, 0.0);
	gsl_matrix_set (m, 1, 0, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 1, 1, -2.0/t);
	gsl_matrix_set (m, 1, 2, 0.0);
	gsl_matrix_set (m, 1, 3, 0.0);
	gsl_matrix_set (m, 2, 0, 0.0);
	gsl_matrix_set (m, 2, 1, 0.0);
	gsl_matrix_set (m, 2, 2, 0.0);
	gsl_matrix_set (m, 2, 3, 1.0);
	gsl_matrix_set (m, 3, 0, -6.0*y[2]*y[0]);
	gsl_matrix_set (m, 3, 1, 0.0);
	gsl_matrix_set (m, 3, 2, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 3, 3, -2.0/t);
	dfdt[0] = 0.0;
	dfdt[1] = 2.0*y[1]/pow(t,2.0);
	dfdt[2] = 0.0;
	dfdt[3] = 2.0*y[3]/pow(t,2.0);
	return GSL_SUCCESS;
	}
	
//gives absolute value of a number
double absolute (const double& amplitude)
	{
	double abs_amplitude;
	if (amplitude > 0)
		{
		abs_amplitude = amplitude;
		}
	else
		{
		abs_amplitude = -amplitude;
		}
	return abs_amplitude;
	}

//simple function to print a vector to file	
void simplePrintVector(const string& printFile, vector<double> vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F.precision(16);
	F << left;
	unsigned int length = vecToPrint.size();
	for (unsigned int j=0; j<length; j++)
		{
		F << setw(25) << j << setw(25) << vecToPrint[j] << endl;
		}
	F.close();
	}
	
//print repi via gnuplot
void gpSimple(const string & readFile) 
	{
	string commandOpenStr = "gnuplot -persistent";
	const char * commandOpen = commandOpenStr.c_str();
	FILE * gnuplotPipe = popen (commandOpen,"w");
	string commandStrPic1 = "set term png size 1600,800";
    string commandStrPic2 = "set output './pics/sphaleron.png'";
	string commandStr1 = "plot \"" + readFile + "\" using 1:2 with linespoints";
	string commandStr2 = "pause -1";
	const char * commandPic1 = commandStrPic1.c_str();
	const char * commandPic2 = commandStrPic2.c_str();
	const char * command1 = commandStr1.c_str();
	const char * command2 = commandStr2.c_str();
	fprintf(gnuplotPipe, "%s \n", commandPic1);
	fprintf(gnuplotPipe, "%s \n", commandPic2);
	fprintf(gnuplotPipe, "%s \n", command1);
	fprintf(gnuplotPipe, "%s \n", command2);
	pclose(gnuplotPipe);
	}
	
struct paramsStruct {double Y_0;};
struct paramsStruct p;
	
//Energy integrand
double EIntegrand (double x, void * parameters)
	{
	if (x<0.0)
		{
		printf ("error in EIntegrand, R<0");
		return 0.0;
		}
	else
		{
		struct paramsStruct * params = (struct paramsStruct *)parameters;
		double Y_0 = (params->Y_0);
		double y_R[4] = { Y_0, 0.0, 1.0, 0.0};
		int status;
		double t = 1.0e-16;
		gsl_odeiv2_system syse = {func, jac, 4, &paramsVoid};
		gsl_odeiv2_driver * de = gsl_odeiv2_driver_alloc_yp_new (&syse,  gsl_odeiv2_step_rk8pd, 1.0e-9, 1.0e-9, 0.0);
		status = gsl_odeiv2_driver_apply (de, &t, x, y_R);
		if (status != GSL_SUCCESS)
			{
			printf ("error in EIntegrand, return value=%d\n", status);
			gsl_odeiv2_driver_free (de);
			return (double)status;
			}
		else
			{
			double to_return = 4.0*pi*pow(x,2.0)*( 0.5*pow(y_R[1],2.0) + 0.5*pow(y_R[0],2.0) - 0.25*pow(y_R[0],4.0) );
			gsl_odeiv2_driver_free (de);
			return to_return;
			}
		}
	}
	
//print sparse matrix to file
void printSpmat (const string & printFile, spMat spmatToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F << left;
	F.precision(16);
	for (int l=0; l<spmatToPrint.outerSize(); ++l)
		{
		for (Eigen::SparseMatrix<double>::InnerIterator it(spmatToPrint,l); it; ++it)
			{
			F << setw(25) << it.row()+1 << setw(25) << it.col()+1 << setw(25) << it.value() << endl;
			}
		}
	F.close();
	}
	
//count non-empty lines of a file
unsigned int countLines(const string & file_to_count)
	{
	ifstream fin;
	fin.open(file_to_count.c_str());
	string line;
	unsigned int counter = 0;
	while(!fin.eof())
		{
		getline(fin,line);
		if(line.empty())
			{
			continue;
			}
		counter++;
		}		
	fin.close();
    return counter;
	}
	
//load simple vector from file
vec loadSimpleVector (const string& loadFile)
	{
	unsigned int fileLength = countLines(loadFile);
	vec outputVec(fileLength);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line;
	unsigned int j=0;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			istringstream ss(line);
			ss >> outputVec(j);
			j++;
			}
		}
	F.close();
	return outputVec;
	}
	
//print three vectors
void printThreeVectors(const string& printFile, vec vec1, vec vec2, vec vec3)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F.precision(20);
	F << left;
	unsigned int length1 = vec1.size();
	unsigned int length2 = vec2.size();
	unsigned int length3 = vec3.size();
	if (length1==length2 && length1==length3)
		{
		for (unsigned int j=0; j<length1; j++)
			{
			F << setw(25) << vec1(j) << setw(25) << vec2(j) << setw(25) << vec3(j) << endl;
			}
		}
	else
		{
		cout << "printThreeVectors error, vectors different lengths: " << length1 << " " << length2 << " " << length3 << endl;
		}
	F.close();
	}
	
//print p via gnuplot, using repi or pi or some such like
void gp(const string & readFile, const string & gnuplotFile) 
	{
	string prefix = "gnuplot -e \"f='";
	string middle = "'\" ";
	string suffix = " -persistent";
	string commandStr = prefix + readFile + middle + gnuplotFile + suffix;
	const char * command = commandStr.c_str();
	FILE * gnuplotPipe = popen (command,"w");
	fprintf(gnuplotPipe, "%s \n", " ");
	pclose(gnuplotPipe);
	}
	
int main()
{
/* first of all just writing something to solve the ode as an initial value plot
only once that works will i add the newton raphson loop to do the shooting
to do:
- define gsl functions for first order odes
- define jacobian functions
- define boundary values
- follow gsl example for implementation

next, have to implement shooting by applying newton's method to F, the difference between the final boundary value and the one shot at
- define function F
- write loop with convergence test
- write newton's method
- print tests for convergence on each iteration

then,
- calculate energy of solution
- find negative mode
- calculate number of particles
	- to do this need to evolve sphaleron + small amount of negative mode in time
	- then model N(k) and E(K) and N_lin and E_lin with time
	- this requires taking something like the fourier transform of the field configuration
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//solving 1d boundary value problem via shooting method
paramsVoid = {};

gsl_odeiv2_system sys = {func, jac, 4, &paramsVoid};

double F = 1.0, dF;
double aim = 0.0;
double closeness = 1.0e-9;
double t0 = 1.0e-16, t1 = 5.0;
const unsigned int N = 2e2;
double h = t1-t0;
h /= (double)N;
vector<double> y0Vec(N+1), y2Vec(N+1);
unsigned int runsCount = 0;

double Y0 = 4.337;//initial guess
//cout << "initial y0: ";
//cin >> Y0;

printf("%16s%16s%16s%16s%16s%16s%16s%16s\n","run","N","y(t1)","yMin","F-aim","Y0Old","Y0New","-F/dF");

while (absolute(F-aim)>closeness)
	{
	runsCount++;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_yp_new (&sys,  gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

	double t = t0, ti;
	double y0[4] = { Y0, 0.0, 1.0, 0.0};
	double y[4];
	memcpy(y, y0, sizeof(y0));
	double yMin = absolute(y[0]-aim);
	y0Vec[0] = y[0];
	y2Vec[0] = y[2];
	int status;
	unsigned int i, iMin = 0;
	
	for (i = 1; i <= N; i++)
		{
		ti =(double)i;
		ti *= h;
		ti += t0;
		status = gsl_odeiv2_driver_apply (d, &t, ti, y);
		if (status != GSL_SUCCESS)
			{
			printf ("error, return value=%d\n", status);
			printf ("i = %3i, t = %3g\n",i,t);
			break;
			}
		if ((y[0]-aim)<yMin && (y[0]-aim)>0.0)
			{
			yMin = y[0]-aim;
			iMin = i;
			}
		y0Vec[i] = y[0];
		y2Vec[i] = y[2];
		//printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
		if ((y[0]-aim)<0.0)
			{
			iMin = i;
			if ((y[0]-aim)<-0.2) { break; }
			}
		}
	if (status != GSL_SUCCESS){break;}
		
	F = y0Vec[iMin]-aim; //as final boundary condition is y(t1)=0.0;
	dF = y2Vec[iMin];
	printf("%16i%16i%16g%16g%16g%16.12g",runsCount,i,y[0],yMin,F-aim,Y0);
	if (absolute(dF)>2.0e-16)
		{
		Y0 += -F/dF;
		printf("%16.12g%16g\n",Y0,-F/dF);
		}
	gsl_odeiv2_driver_free (d);
	if (i==(N+1))
		{
		F = y[0];
		}
	}

//printing solution
string filename = "./data/sphaleron.dat";
simplePrintVector(filename,y0Vec);
printf("\n%16s%20s%20s\n\n","Solution printed: ",filename.c_str()," ./pics/sphaleron.png");
gpSimple(filename);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//finding E
double E, Eerror;
p  = {Y0};
gsl_function E_integrand;
E_integrand.function = &EIntegrand;
E_integrand.params = &p;
gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
gsl_integration_qag(&E_integrand, t0, t1, 1.0e-10, 1.0e-9, 1e4, 4, w, &E, &Eerror);
gsl_integration_workspace_free(w);
if (Eerror>1.0e-8) { cout << "E error = " << Eerror << endl;}
else { cout << "E = " << E << endl;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//computing and printing linear fluctuation operator as sparse matrix
//in order to calulate negative eigenvalue and eigenvector
spMat D1(N+1,N+1), D2(N+1,N+1);
D1.setZero();
D2.setZero();
double t = t0;
for (unsigned int j=0; j<(N+1); j++)
	{
	if (j==0)
		{
		D1.insert(j,j) = pow(t,2.0)/h + pow(t,2.0)*h*( 1.0 - 3.0*pow(y0Vec[j],2.0) );
		D1.insert(j,j+1) = -pow(t,2.0)/h;
		D2.insert(j,j) = 1.0/pow(h,2.0) + 2.0/h/t + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D2.insert(j,j+1) = -1.0/pow(h,2.0) - 2.0/h/t;
		/*D1.insert(j,j) = -1.0; //derivative of phi is zero at r=0
		D1.insert(j,j+1) = 1.0;
		D2.insert(j,j) = -1.0;
		D2.insert(j,j+1) = 1.0;*/
		}
	else if (j==N)
		{
		D1.insert(j,j) = pow(t-h,2.0)/h + pow(t,2.0)*h*( 1.0 - 3.0*pow(y0Vec[j],2.0) );
		D1.insert(j,j-1) = -pow(t-h,2.0)/h;
		D2.insert(j,j) = 1.0/pow(h,2.0) + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D2.insert(j,j-1) = -1.0/pow(h,2.0);
		/*D1.insert(j,j) = 1.0; //phi goes to zero as r->infty
		D2.insert(j,j) = 1.0;*/
		}
	else
		{
		D1.insert(j,j) = pow(t,2.0)/h + pow(t-h,2.0)/h + pow(t,2.0)*h*( 1.0 - 3.0*pow(y0Vec[j],2.0) );
		D1.insert(j,j+1) = -pow(t,2.0)/h;
		D1.insert(j,j-1) = -pow(t-h,2.0)/h;
		D2.insert(j,j) = 2.0/pow(h,2.0) + 2.0/h/t + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D2.insert(j,j+1) = -1.0/pow(h,2.0) - 2.0/h/t;
		D2.insert(j,j-1) = -1.0/pow(h,2.0);
		}
	t += h;
	}
D1.makeCompressed();
D2.makeCompressed();
printSpmat("./data/D1.dat",D1);
printSpmat("./data/D2.dat",D2);
printf("\nMatrices printed: ./data/D1.dat ./data/D2.dat\n\n");
printf("From Matlab, omega^2_- = -15.3060\n\n");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//propagating sphaleron plus a small amount of negative mode forward in time
//in order to calculate E_lin and N_lin

unsigned int Nt = 1e3;
double T = 10.0, amp = 1.0e-3;
double dt = T/Nt; //equals Dt
vec phi((Nt+1)*(N+1)), vel((Nt+1)*(N+1)), acc((Nt+1)*(N+1)), eigVec, xVec((Nt+1)*(N+1)), tVec((Nt+1)*(N+1));

//getting eigVec from file
eigVec = loadSimpleVector("data/sphaleronEigVec.dat");

//defining phi on initial time slice
for (unsigned int j=0;j<(N+1);j++)
	{
	unsigned int l = j*(Nt+1);
	tVec(l) = 0.0;
	xVec(l) = t0+j*h;
	phi(l) = y0Vec[j] + amp*eigVec(j);
	}

//intitialize velocity
for (unsigned int j=0; j<(N+1); j++)
	{
	unsigned int l = j*(Nt+1);
    vel(l) = 0.0; //not sure about this initial condition
	}


//initialize acc using phi and expression from equation of motion
acc(0) = (phi((Nt+1)) - phi(0))/pow(h,2.0) - phi(0) + pow(phi(0),3.0);
unsigned int l = N*(Nt+1);
acc(l) = (phi(l-(Nt+1))*t1/(t1-h) - phi(l)*t1/(t1-h))/pow(h,2.0) - phi(l) + pow(phi(l),3.0);
for (unsigned int j=1; j<N; j++)
	{
	l = j*(Nt+1);
	double r = t0 + h*j;
    acc(l) = (phi(l+(Nt+1)) + phi(l-(Nt+1))*r/(r-h) - phi(l)*(1.0+r/(r-h)))/pow(h,2.0) - phi(l) + pow(phi(l),3.0);
	}
	
//A4.5 starting the energy and particle number off
//vec linErg(Nt+1); linErg = Eigen::VectorXd::Zero(Nt+1);
//vec linNum(Nt+1); linNum = Eigen::VectorXd::Zero(Nt+1);

//A7. run loop
for (unsigned int u=1; u<(Nt+1); u++)
	{
    for (unsigned int x=0; x<(N+1); x++)
    	{
        unsigned int m = u+x*(Nt+1);
        tVec(l) = u*dt;
		xVec(l) = t0+x*h;
        vel(m) = vel(m-1) + dt*acc(m-1);
        phi(m) = phi(m-1) + dt*vel(m);
    	}
    acc(u) = (phi(u+(Nt+1)) - phi(u))/pow(h,2.0) - phi(u) + pow(phi(u),3.0);
    unsigned int v = u+N*(Nt+1);
    acc(v) = (phi(v-(Nt+1))*t1/(t1-h) - phi(v)*t1/(t1-h))/pow(h,2.0) - phi(v) + pow(phi(v),3.0);
    for (unsigned int x=1; x<N; x++)
    	{
        unsigned int m = u+x*(Nt+1);
        double r = t0 + x*h;
        acc(m) = (phi(m+(Nt+1)) + phi(m-(Nt+1))*r/(r-h) - phi(m)*(1.0+r/(r-h)))/pow(h,2.0) - phi(m) + pow(phi(m),3.0);
    	}
	}
	
printThreeVectors("data/sphaleronEvolution.dat",tVec,xVec,phi);
gp("data/sphaleronEvolution.dat","sphaleron.gp");
printf("Time evolution printed: data/sphaleronEvolution.dat pics/sphaleronEvolution.png\n");

return 0;
}
