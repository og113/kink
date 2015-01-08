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

struct force_params{double x0; double x1; unsigned int Nx; double sigma;}; 
struct force_params paramsForce;

#define pi 3.14159265359

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//main ode and pde functions

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
	gsl_matrix_set_all (m, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 1, 0, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 1, 1, -2.0/t);
	gsl_matrix_set (m, 2, 3, 1.0);
	gsl_matrix_set (m, 3, 0, -6.0*y[2]*y[0]);
	gsl_matrix_set (m, 3, 2, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 3, 3, -2.0/t);
	dfdt[0] = 0.0;
	dfdt[1] = 2.0*y[1]/pow(t,2.0);
	dfdt[2] = 0.0;
	dfdt[3] = 2.0*y[3]/pow(t,2.0);
	return GSL_SUCCESS;
	}
	
//a function to give 'force' for time evolution
int funcPDE (double t, const double y[], double f[], void *params)
	{
	struct force_params * parameters = (struct force_params *)params;
	double x0 = (parameters->x0);
	double x1 = (parameters->x1);
	unsigned int Nx = (parameters->Nx);
	double sigma = (parameters->sigma);
	double dx = (x1-x0), x;
	dx /= (double)Nx;
	unsigned int yLength = (unsigned int)(sizeof(y)/sizeof(*y));
	if (yLength!=2*Nx){printf("funcPDE error: yLength = %i",yLength);}
	f[0] = y[1];
	f[1] = 2.0*(y[2]-y[0])/pow(dx,2.0) - y[0] + pow(y[0],3.0);
	f[1] *= sigma;
	for (unsigned int j=1; j<(Nx-1); j++)
		{
		x = x0 + (double)j*dx;
		f[2*j] = y[2*j+1];
		f[2*j+1] = (y[2*(j+1)] + y[2*(j-1)]-2.0*y[2*j])/pow(dx,2.0) + (y[2*(j+1)] - y[2*(j-1)])/x/dx - y[2*j] + pow(y[2*j],3.0);
		f[2*j+1] *= sigma;
		}
	x = x0 + (double)(Nx-1.0)*dx;
	f[2*(Nx-1)] = y[2*(Nx-1)+1];
	f[2*(Nx-1)+1] = (y[2*(Nx-2)]-2.0*y[2*(Nx-1)])/pow(dx,2.0) + (- y[2*(Nx-2)])/x/dx - y[2*(Nx-1)] + pow(y[2*(Nx-1)],3.0);
	f[2*(Nx-1)+1] *= sigma;
	return GSL_SUCCESS;
	}
	
//a function to give 'force' for time evolution	
int jacPDE (double t, const double y[], double *dfdy, double dfdt[], void *params)
	{
	struct force_params * parameters = (struct force_params *)params;
	double x0 = (parameters->x0);
	double x1 = (parameters->x1);
	unsigned int Nx = (parameters->Nx);
	double sigma = (parameters->sigma);
	double dx = (x1-x0), x;
	dx /= (double)Nx;
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2*Nx, 2*Nx);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set_all (m, 0.0);
	for (unsigned int j=0; j<Nx; j++)
		{
		if (j==0)
			{
			gsl_matrix_set(m,2*j,2*j+1, 1.0 );
			gsl_matrix_set(m,2*j+1,2*j, sigma*( -2.0/pow(dx,2.0) - 1.0 + 3.0*pow(y[2*j],2.0) ) );
			gsl_matrix_set(m,2*j+1,2*(j+1), sigma*( 2.0/pow(dx,2.0) ) );
			dfdt[2*j] = 0.0;
			dfdt[2*j+1] = 0.0;
			}
		else if (j==(Nx-1))
			{
			x = x0 + (double)(Nx-1.0)*dx;
			gsl_matrix_set(m,2*j,2*j+1, 1.0 );
			gsl_matrix_set(m,2*j+1,2*j, sigma*( -2.0/pow(dx,2.0) - 1.0 + 3.0*pow(y[2*j],2.0) ) );
			gsl_matrix_set(m,2*j+1,2*(j-1), sigma*( 1.0/pow(dx,2.0) - 1.0/x/dx ) );
			dfdt[2*j] = 0.0;
			dfdt[2*j+1] = 0.0;
			}
		else
			{
			x = x0 + (double)j*dx;
			gsl_matrix_set(m,2*j,2*j+1, 1.0 );
			gsl_matrix_set(m,2*j+1,2*j, sigma*( -2.0/pow(dx,2.0) - 1.0 + 3.0*pow(y[2*j],2.0) ) );
			gsl_matrix_set(m,2*j+1,2*(j+1), sigma*( 1.0/pow(dx,2.0) + 1.0/x/dx ) );
			gsl_matrix_set(m,2*j+1,2*(j-1), sigma*( 1.0/pow(dx,2.0) - 1.0/x/dx ) );
			dfdt[2*j] = 0.0;
			dfdt[2*j+1] = 0.0;
			}
		}
	return GSL_SUCCESS;
	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
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

//some simple vector functions
vector<double> vecEquate(const vector<double> & vecB)
	{
	unsigned int length = vecB.size();
	vector<double> vecA(length);
	for (unsigned int j=0; j<length; j++)
		{
		vecA[j] = vecB[j];
		}
	return vecA;
	}

vector<double> vecSum(const vector<double> & vecB, const vector<double> & vecA)
	{
	unsigned int lengthA = vecA.size();
	unsigned int lengthB = vecB.size();
	if (lengthA!=lengthB)
		{
		printf("vecSum error: vector sizes do not agree: lengthA = %i, lengthB = %i\n",lengthA,lengthB);
		return vecA;
		}
	vector<double> vecC(lengthA);
	for (unsigned int j=0; j<lengthA; j++)
		{
		vecC[j] = vecA[j] + vecB[j];
		}
	return vecC;
	}
	
vector<double> vecMultiply(const vector<double> & vecA, const double alpha)
	{
	unsigned int lengthA = vecA.size();
	vector<double> vecB(lengthA);
	for (unsigned int j=0; j<lengthA; j++)
		{
		vecB[j] = alpha*vecA[j];
		}
	return vecB;
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
		F << setw(25) << j << setw(25) << vecToPrint[j] << endl;
		}
	F.close();
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
		printf ("error in EIntegrand, R<0\n");
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
	
//print matrix to file
void printMat (const string & printFile, mat matToPrint)
	{
	fstream F;
	F.open((printFile).c_str(), ios::out);
	F << left;
	F.precision(16);
	for (unsigned int l=0; l<matToPrint.rows(); l++)
		{
		for (unsigned int m=0; m<matToPrint.cols(); m++)
			{
			F << setw(25) << matToPrint(l,m) << endl;
			}
		}
	F.close();
	}
	
//load matrix from file
mat loadMat (const string & loadFile, const unsigned int& rows, const unsigned int & cols)
	{
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	mat matOut(rows,cols);
	for (unsigned int l=0; l<rows; l++)
		{
		for (unsigned int m=0; m<cols; m++)
			{
			F >> matOut(l,m);
			}
		}
	F.close();
	return matOut;
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
			ss >> outputVec[j];
			j++;
			}
		}
	F.close();
	return outputVec;
	}
	
//load simple vector from file
vec loadSimpleVectorColumn (const string& loadFile, const unsigned int column)
	{
	unsigned int fileLength = countLines(loadFile);
	vec outputVec(fileLength);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line, temp;
	unsigned int j=0;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			istringstream ss(line);
			for (unsigned int l=0; l<column; l++)
				{
				ss >> temp;
				}
			ss >> outputVec[j];
			j++;
			}
		}
	F.close();
	return outputVec;
	}
	
//load simple vector from file
vector<double> loadSimpleVector2 (const string& loadFile)
	{
	unsigned int fileLength = countLines(loadFile);
	vector<double> outputVec(fileLength);
	fstream F;
	F.open((loadFile).c_str(), ios::in);
	string line;
	unsigned int j=0;
	while (getline(F, line))
		{
		if (!line.empty())
			{
			istringstream ss(line);
			ss >> outputVec[j];
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
	
//same function but for vector<double>	
vector<double> interpolate2(vector<double> vec_old, const unsigned int & Nt_old, const unsigned int & N_old, const unsigned int & Nt_new, const unsigned int & N_new)
	{
	unsigned int old_size = vec_old.size();
	if (old_size<N_old*Nt_old) {cout << "interpolate error, vec_old.size() = " << old_size << " , N_old*Nt_old = " << N_old*Nt_old << endl;}
	vector<double> vec_new (N_new*Nt_new);
	
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
			vec_new[l] = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old[pos] \
							+ (1.0-rem_t_old)*rem_x_old*vec_old[pos+Nt_old] \
							+ rem_t_old*(1.0-rem_x_old)*vec_old[pos+1] \
							+ rem_t_old*rem_x_old*vec_old[pos+Nt_old+1];
			}
		else if (x_old<(N_old-1))
			{
			vec_new[l] = (1.0-rem_x_old)*vec_old[pos] + rem_x_old*vec_old[pos+Nt_old];
			}
		else
			{
			vec_new[l] = vec_old[pos];
			}
		}
	return vec_new;
	}
	
//same function but 1d	
vec interpolate1d(vec vec_old, const unsigned int & N_old, const unsigned int & N_new)
	{
	unsigned int old_size = vec_old.size();
	if (old_size<N_old) {cout << "interpolate error, vec_old.size() = " << old_size << " , N_old = " << N_old << endl;}
	vec vec_new (N_new);
	
	unsigned int x_old;
	double exact_x_old, rem_x_old;
	
	for (unsigned int l=0;l<N_new;l++)
		{
		exact_x_old = l*(N_old-1.0)/(N_new-1.0);
		x_old = (unsigned int)exact_x_old;
		rem_x_old = exact_x_old - (double)(x_old);
		if  (x_old<(N_old-1))
			{
			vec_new[l] = (1.0-rem_x_old)*vec_old[x_old] \
							+ rem_x_old*vec_old[x_old+1];
			}
		else
			{
			vec_new[l] = vec_old[x_old];
			}
		}
	return vec_new;
	}
	
//a function to make a vector defined in 2d equal to one defined in 1d on a given timeslice
void equateSlice(vec vec2d, vec vec1d, const unsigned int & Ntime, const unsigned int & Nspace, const unsigned int & tSlice)
	{
	if (vec2d.size()!=Ntime*Nspace || vec1d.size()!=Nspace)
		{
		printf("equateSlice error: vec2d.length() = %i, vec1d.length() = %i\nNtime = %i, Nspace = %i\n",(int)(vec2d.size()),(int)(vec1d.size()),Ntime,Nspace);
		}
	for (unsigned int l=0;l<Nspace;l++)
		{
		unsigned int m = tSlice + Ntime*l;
		vec2d(m) = vec1d(l);
		}
	}
	
//function to reverse time direction of a vector
vec reverseTime(vec inVec, const unsigned int & Nt, const unsigned int & Nx)
	{
	vec outVec(Nx*Nt);
	for (unsigned int j=0; j<Nx*Nt; j++)
		{
		unsigned int x = intCoord(j,1,Nt), t = intCoord(j,0,Nt);
		unsigned int m = t + x*Nt, n = (Nt-1-t) + x*Nt;
		outVec[n] = inVec[m];
		}
	return outVec;
	}
	
//getting the date and time
const string currentDateTime()
	{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%y%m%d%H%M%S", &tstruct);
    return buf;
	}
