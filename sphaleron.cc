//parameters and functions for pi.cc
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
void simplePrintVector(const string& printFile, vec vecToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F.precision(16);
	F << left;
	unsigned int length = vecToPrint.size();
	for (unsigned int j=0; j<length; j++)
		{
		F << setw(25) << j << setw(25) << vecToPrint(j) << endl;
		}
	F.close();
	}
	
//print repi via gnuplot
void gpSimple(const string & readFile, const string & outFile) 
	{
	string commandOpenStr = "gnuplot -persistent";
	const char * commandOpen = commandOpenStr.c_str();
	FILE * gnuplotPipe = popen (commandOpen,"w");
	string commandStrPic1 = "set term png size 1600,800";
    string commandStrPic2 = "set output '" + outFile + "'";
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
	
//function which find integer coordinates in 2d
unsigned int intCoord(const unsigned int& locNum, const int& direction, const unsigned int& xNt)
	{
	unsigned int XintCoord;
	unsigned int x = floor(locNum/xNt);
	if (direction==1)
		{
		XintCoord = x;
		}
	else
		{
		XintCoord = locNum - x*xNt;
		}
	return XintCoord;	
	}
	
//interpolate vector function
vec interpolate(vec vec_old, const unsigned int & Nt_old, const unsigned int & N_old, const unsigned int & Nt_new, const unsigned int & N_new)
	{
	unsigned int old_size = vec_old.size();
	if (old_size<N_old*Nt_old) {cout << "interpolate error, vec_old.size() = " << old_size << " , N_old*Nt_old = " << N_old*Nt_old << endl;}
	vec vec_new (N_new*Nt_new);
	
	unsigned int x_new, t_new, x_old, t_old;
	double exact_x_old, exact_t_old, rem_x_old, rem_t_old;
	unsigned int pos;	
	
	for (unsigned int l=0;l<N_new*Nt_new;l++)
		{
		t_new = intCoord(l,0,Nt_new);
		x_new = intCoord(l,1,Nt_new);
		exact_t_old = t_new*(Nt_old-1.0)/(Nt_new-1.0);
		exact_x_old = x_new*(N_old-1.0)/(N_new-1.0);
		t_old = (unsigned int)exact_t_old;
		x_old = (unsigned int)exact_x_old;
		rem_t_old = exact_t_old;
		rem_t_old -= (double)(t_old);
		rem_x_old = exact_x_old;
		rem_x_old -= (double)(x_old);
		pos = t_old + Nt_old*x_old;
		if  (t_old<(Nt_old-1) && x_old<(N_old-1))
			{
			vec_new(l) = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old(pos) \
							+ (1.0-rem_t_old)*rem_x_old*vec_old(pos+Nt_old) \
							+ rem_t_old*(1.0-rem_x_old)*vec_old(pos+1) \
							+ rem_t_old*rem_x_old*vec_old(pos+Nt_old+1);
			}
		else if (x_old<(N_old-1))
			{
			vec_new(l) = (1.0-rem_x_old)*vec_old(pos) + rem_x_old*vec_old(pos+Nt_old);
			}
		else
			{
			vec_new(l) = vec_old(pos);
			}
		}
	return vec_new;
	}
	
int main()
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//solving 1d boundary value problem via shooting method
paramsVoid = {};

gsl_odeiv2_system sys = {func, jac, 4, &paramsVoid};

double F = 1.0, dF;
double aim = 0.0;
double closeness = 1.0e-9;
double r0 = 1.0e-16, r1 = 5.0;
const unsigned int N = 1e4;
double dr = r1-r0;
dr /= (double)N;
vec y0Vec(N+1), y2Vec(N+1);
unsigned int runsCount = 0;

double Y0 = 4.337;//initial guess
//cout << "initial y0: ";
//cin >> Y0;

printf("%16s%16s%16s%16s%16s%16s%16s%16s\n","run","N","y(r1)","yMin","F-aim","Y0Old","Y0New","-F/dF");

while (absolute(F-aim)>closeness)
	{
	runsCount++;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_yp_new (&sys,  gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

	double r = r0, ri;
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
		ri =(double)i;
		ri *= dr;
		ri += r0;
		status = gsl_odeiv2_driver_apply (d, &r, ri, y);
		if (status != GSL_SUCCESS)
			{
			printf ("error, return value=%d\n", status);
			printf ("i = %3i, r = %3g\n",i,r);
			break;
			}
		if ((y[0]-aim)<yMin && (y[0]-aim)>0.0)
			{
			yMin = y[0]-aim;
			iMin = i;
			}
		y0Vec[i] = y[0];
		y2Vec[i] = y[2];
		//printf ("%.5e %.5e %.5e\n", r, y[0], y[1]);
		if ((y[0]-aim)<0.0)
			{
			iMin = i;
			if ((y[0]-aim)<-0.2) { break; }
			}
		}
	if (status != GSL_SUCCESS){break;}
		
	F = y0Vec[iMin]-aim; //as final boundary condition is y(r1)=0.0;
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
string filename = "./data/sphaleron.dat", picname = "./pics/sphaleron.png";
simplePrintVector(filename,y0Vec);
printf("\n%16s%20s%20s\n\n","Solution printed: ",filename.c_str(),picname.c_str());
gpSimple(filename,picname);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//finding E
double E, Eerror;
p  = {Y0};
gsl_function E_integrand;
E_integrand.function = &EIntegrand;
E_integrand.params = &p;
gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
gsl_integration_qag(&E_integrand, r0, r1, 1.0e-10, 1.0e-9, 1e4, 4, w, &E, &Eerror);
gsl_integration_workspace_free(w);
if (Eerror>1.0e-8) { cout << "E error = " << Eerror << endl;}
else { cout << "E = " << E << endl;}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//computing and printing linear fluctuation operator as sparse matrix
//in order to calulate negative eigenvalue and eigenvector
spMat D1(N+1,N+1), D2(N+1,N+1);
D1.setZero();
D2.setZero();
double r = r0;
for (unsigned int j=0; j<(N+1); j++)
	{
	if (j==0)
		{
		D1.insert(j,j) = 1.0/pow(dr,2.0) + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D1.insert(j,j+1) = -1.0/pow(dr,2.0);
		D2.insert(j,j) = 1.0/pow(dr,2.0) + 2.0/dr/r + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D2.insert(j,j+1) = -1.0/pow(dr,2.0) - 2.0/dr/r;
		/*D1.insert(j,j) = -1.0; //derivative of phi is zero at r=0
		D1.insert(j,j+1) = 1.0;
		D2.insert(j,j) = -1.0;
		D2.insert(j,j+1) = 1.0;*/
		}
	else if (j==N)
		{
		D1.insert(j,j) = 1.0/pow(dr,2.0) - 2.0*(1.0-dr/r)/r/dr + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D1.insert(j,j-1) = -1.0/pow(dr,2.0) + 2.0*(1.0-dr/r)/r/dr;
		D2.insert(j,j) = 1.0/pow(dr,2.0) + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D2.insert(j,j-1) = -1.0/pow(dr,2.0);
		/*D1.insert(j,j) = 1.0; //phi goes to zero as r->infty
		D2.insert(j,j) = 1.0;*/
		}
	else
		{
		D1.insert(j,j) = 2.0/pow(dr,2.0) - 2.0*(1.0-dr/r)/r/dr + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D1.insert(j,j+1) = -1.0/pow(dr,2.0);
		D1.insert(j,j-1) = -1.0/pow(dr,2.0) + 2.0*(1.0-dr/r)/r/dr;
		D2.insert(j,j) = 2.0/pow(dr,2.0) + 2.0/dr/r + 1.0 - 3.0*pow(y0Vec[j],2.0);
		D2.insert(j,j+1) = -1.0/pow(dr,2.0) - 2.0/dr/r;
		D2.insert(j,j-1) = -1.0/pow(dr,2.0);
		}
	r += dr;
	}
D1.makeCompressed();
D2.makeCompressed();
printSpmat("./data/D1.dat",D1);
printSpmat("./data/D2.dat",D2);
printf("\nMatrices printed: ./data/D1.dat ./data/D2.dat\n\n");
printf("From Matlab, omega^2_- = -15.3060\n\n");
return 0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//assembling material to calculate occupation of modes

//h_ij matrix
mat h(N+1,N+1);
for (unsigned int j=0;j<(N+1);j++)
	{
	if (j==0)
		{
		h(j,j) = 1.0 + 1.0/pow(dr,2.0);
		h(j,j+1) = -1.0/pow(dr,2.0);
		}
	else if (j==N)
		{
		h(j,j) = 1.0/pow(dr,2.0);
		h(j,j-1) = -1.0/pow(dr,2.0);
		}
	else
		{
		h(j,j) = 1.0 + 2.0/pow(dr,2.0);
		h(j,j+1) = -1.0/pow(dr,2.0);
		h(j,j-1) = -1.0/pow(dr,2.0);
		}
	}
mat omega(N,N); 	omega = Eigen::MatrixXd::Zero(N,N);
mat Eomega(N,N); 	Eomega = Eigen::MatrixXd::Zero(N,N);
vec eigenValues(N);
mat eigenVectors(N,N); //eigenvectors correspond to columns of this matrix
Eigen::SelfAdjointEigenSolver<mat> eigensolver(h);
if (eigensolver.info() != Eigen::Success)
	{
	cout << "h eigensolver failed" << endl;
	}
else
	{
	eigenValues = eigensolver.eigenvalues();
	eigenVectors = eigensolver.eigenvectors(); //automatically normalised to have unit norm
	}
	
for (unsigned int j=0; j<N; j++)
	{
	for (unsigned int k=0; k<N; k++)
		{
		double dj_sqrt = (r0 + j*dr)*pow(dr,0.5);
		double dk_sqrt = (r0 + k*dr)*pow(dr,0.5);
		for (unsigned int l=0; l<N; l++)
			{
			omega(j,k) += dj_sqrt*dk_sqrt*pow(eigenValues(l),0.5)*eigenVectors(j,l)*eigenVectors(k,l);
			Eomega(j,k) += dj_sqrt*dk_sqrt*eigenValues(l)*eigenVectors(j,l)*eigenVectors(k,l);
			}
		}
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//propagating sphaleron plus a small amount of negative mode forward in time
//in order to calculate E_lin and N_lin

unsigned int Nt = N*3;
double T = 16.0, amp = -5.0e-3;
if ((T-4.0)>1.1*(r1-r0)) { cout << "R is too small compared to T. R = " << r1-r0 << ", T = " << T << endl;}
//cout << "amp of negative mode: ";
//cin >> amp;
double dt = T/Nt; //equals Dt
if (dt>0.5*dr)
	{
	cout << "dt too large. dt = "<< dt << ", dr = " << dr << endl;
	return 0;
	}
vec phi((Nt+1)*(N+1)), vel((Nt+1)*(N+1)), acc((Nt+1)*(N+1)), eigVec, linNum(Nt+1), linErg(Nt+1);
linErg = Eigen::VectorXd::Zero(Nt);
linNum = Eigen::VectorXd::Zero(Nt);

//getting eigVec from file
eigVec = loadSimpleVector("data/sphaleronEigVec.dat");

//initial condition 1)
//defining phi on initial time slice
for (unsigned int j=0;j<(N+1);j++)
	{
	unsigned int l = j*(Nt+1);
	phi(l) = y0Vec[j] + amp*eigVec(j);
	}
	
//initialising linErg and linNum	
for (unsigned int x=0;x<N;x++)
	{
	for (unsigned int k=0;k<N;k++)
		{
		unsigned int j = x*Nt;
		unsigned int l = k*Nt;
		linErg(0) += Eomega(x,k)*phi(l)*phi(j);
		linNum(0) += omega(x,k)*phi(l)*phi(j);
		}
	}

//initial condition 2)
//intitialize velocity, dphi/dt(x,0) = 0;
for (unsigned int j=0; j<(N+1); j++)
	{
	unsigned int l = j*(Nt+1);
    vel(l) = 0.0;
	}

//boundary condition 1)
//phi=0.0 at r=L
for (unsigned int j=0;j<(Nt+1);j++)
	{
	unsigned int l = N*(Nt+1) + j;
	phi(l) = 0.0;
	}

//boundary condition 2)
//initialize acc using phi and expression from equation of motion
/*the unusual d^2phi/dr^2 term and the absence of the first derivative term
are due to the boundary condition 2) which, to second order, is phi(t,-dr) = phi(t,dr)*/
//acc(0) = 2.0*(phi(Nt+1) - phi(0))/pow(dr,2.0) - phi(0) + pow(phi(0),3.0);
acc(0) *= 0.5; //as initial time slice, generated from taylor expansion and equation of motion
for (unsigned int j=1; j<N; j++)
	{
	unsigned int l = j*(Nt+1);
	double r = r0 + dr*j;
    //acc(l) = (phi(l+(Nt+1)) + phi(l-(Nt+1)) - 2.0*phi(l))/pow(dr,2.0) + (phi(l+(Nt+1))-phi(l-(Nt+1)))/r/dr - phi(l) + pow(phi(l),3.0);
    acc(l) = (phi(l+(Nt+1)) + phi(l-(Nt+1)) - 2.0*phi(l))/pow(dr,2.0) + 2.0*(r-dr)*(phi(l)-phi(l-(Nt+1)))/pow(r,2.0)/dr - phi(l) + pow(phi(l),3.0);
    acc(l) *= 0.5;
	}
	
//A4.5 starting the energy and particle number off
//vec linErg(Nt+1); linErg = Eigen::VectorXd::Zero(Nt+1);
//vec linNum(Nt+1); linNum = Eigen::VectorXd::Zero(Nt+1);

//A7. run loop
for (unsigned int u=1; u<(Nt+1); u++)
	{
    for (unsigned int x=0; x<N; x++) //don't loop over last x position as fixed by boundary condition 1)
    	{
        unsigned int m = u+x*(Nt+1);
        vel(m) = vel(m-1) + dt*acc(m-1);
        phi(m) = phi(m-1) + dt*vel(m);
    	}
    //acc(u) = 2.0*(phi(u+Nt+1) - phi(u))/pow(dr,2.0) - phi(u) + pow(phi(u),3.0);
    acc(u) = 2.0*(phi(u+Nt+1) - phi(u))/pow(dr,2.0) - phi(u) + pow(phi(u),3.0);
    for (unsigned int x=1; x<N; x++)
    	{
        unsigned int m = u+x*(Nt+1);
        double r = r0 + x*dr;
        //acc(m) = (phi(m+(Nt+1)) + phi(m-(Nt+1)) - 2.0*phi(m))/pow(dr,2.0) + (phi(m+(Nt+1))-phi(m-(Nt+1)))/r/dr - phi(m) + pow(phi(m),3.0);
        acc(m) = (phi(m+(Nt+1)) + phi(m-(Nt+1)) - 2.0*phi(m))/pow(dr,2.0) + 2.0*(r-dr)*(phi(m)-phi(m-(Nt+1)))/pow(r,2.0)/dr - phi(m) + pow(phi(m),3.0);
        for (unsigned int k=0;k<N;k++)
			{
			unsigned int j = u+k*Nt;
			linErg(u) += Eomega(x,k)*phi(m)*phi(j);
			linNum(u) += omega(x,k)*phi(m)*phi(j);
			}
    	}
	}
	
unsigned int N_print = 100, Nt_print = 100;
vec tVec(Nt_print*N_print), rVec(Nt_print*N_print);
double dtPrint = T/(Nt_print-1.0);
double dxPrint = (r1-r0)/(N_print-1.0);
for (unsigned int t=0;t<Nt_print;t++)
	{
	for (unsigned int r=0; r<N_print; r++)
		{
		unsigned int j= t + r*Nt_print;
		tVec(j) = t*dtPrint;
		rVec(j) = r0 + r*dxPrint;
		}
	}
	
vec phiToPrint;
phiToPrint = interpolate(phi,Nt+1,N+1,Nt_print,N_print);
simplePrintVector("data/sphaleronLinNum.dat",linNum);
simplePrintVector("data/sphaleronLinErg.dat",linErg);
gpSimple("data/sphaleronLinNum.dat","pics/sphaleronLinNum.png");
gpSimple("data/sphaleronLinErg.dat","pics/sphaleronLinErg.png");
printThreeVectors("data/sphaleronEvolution.dat",tVec,rVec,phiToPrint);
gp("data/sphaleronEvolution.dat","sphaleron.gp");
printf("Time evolution printed: data/sphaleronEvolution.dat pics/sphaleronEvolution.png\n");
printf("                        data/sphaleronLinNum.dat pics/sphaleronLinNum.png\n");
printf("                        data/sphaleronLinErg.dat pics/sphaleronLinErg.png\n");

return 0;
}
