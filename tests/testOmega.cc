/*
	program to test omega in various different formulations
*/

#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>

using namespace std;

typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;

#define pi 3.14159265359

// pair for sorting eigenvalues
typedef std::pair<double,size_t> mypair;
bool comparator ( const mypair& l, const mypair& r)
   { return l.first < r.first; }

//to convert string to number, usage is (number type) x = StringToNumber<number type>(str);
template <typename T>
T stringToNumber ( const string &Text ){
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

//gives absolute measure of difference between two numbers
double absDiff (const double& numA, const double& numB) {
	return 2.0*abs(numA-numB)/sqrt(numA*numA+numB*numB);
}

//gives absolute measure of difference between two vectors
double absDiff (const vec& vecA, const vec& vecB) {
	vec diff = vecA-vecB;
	double normDiff = diff.norm();
	double squaredNormA = vecA.squaredNorm();
	double squaredNormB = vecB.squaredNorm();
	return 2.0*normDiff/sqrt(squaredNormA+squaredNormB);
}

//h the matrix from dl[7]
mat hFn(const unsigned int & xN, const double & xa, const double & xmass2)
	{
	mat xh(xN,xN);	xh = Eigen::MatrixXd::Zero(xN,xN);
	double diag = xmass2 + 2.0/pow(xa,2.0);
	double offDiag1 = -1.0/pow(xa,2.0);
	double offDiag2 = -pow(2.0,0.5)/pow(xa,2.0);	
	for (unsigned int l=0; l<xN; l++)
		{
		if (l==0)
			{
			xh(l,l) = 1.0;
			//xh(l,l) = diag;
			//xh(l,l+1) = offDiag2;			
			}
		else if (l==(xN-1))
			{
			xh(l,l) = 1.0;
			//xh(l,l) = diag;
			//xh(l,l-1) = offDiag2;
			}
		else
			{
			xh(l,l) = diag;
			xh(l,l+1) = ((l+1)==(xN-1)? offDiag2: offDiag1);
			xh(l,l-1) = ((l-1)==0?		offDiag2: offDiag1);
			}
		}
	return xh;
	}

int main(int argc, char** argv) {

/*---------------------------------------------------------------------------------------------
	main parameters
---------------------------------------------------------------------------------------------*/
unsigned int N = 10;
double L = 1.0;

/*---------------------------------------------------------------------------------------------
	getting inputss
---------------------------------------------------------------------------------------------*/

if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("L")==0) L = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("N")==0) N = atoi(argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

mat omega1(N,N), omega2(N,N), Eomega1(N,N), Eomega2(N,N);
double dr = L/(N-1.0);

printf("%12s%12s%12s\n","N","L","dr");
printf("%12i%12g%12g\n\n",N,L,dr);

/*---------------------------------------------------------------------------------------------
	1) expression based on numerical calculation of eigenvectors
---------------------------------------------------------------------------------------------*/
mat h(N,N);
h = hFn(N,dr,1.0);
//cout << "h = " << endl << h << endl << endl;
omega1 = Eigen::MatrixXd::Zero(N,N);
Eomega1 = Eigen::MatrixXd::Zero(N,N);
vec eigenValues1(N);
mat eigenVectors1(N,N);
Eigen::VectorXcd eigenValues1Complex(N);
Eigen::MatrixXcd eigenVectors1Copy(N,N), eigenVectors1Complex(N,N); //eigenvectors correspond to columns of this matrix
Eigen::EigenSolver<mat> eigensolver(h);
if (eigensolver.info() != Eigen::Success) cerr << "h eigensolver failed" << endl;
else {
	eigenValues1Complex = eigensolver.eigenvalues();
	eigenVectors1Copy = eigensolver.eigenvectors(); //automatically normalised to have unit norm
}
//cout << eigenValues1Complex << endl << endl;
vector<mypair> eigenvaluesPair(N);
for (unsigned int j=0; j<N; j++) {
	eigenvaluesPair[j].first = real(eigenValues1Complex[j]);
	eigenvaluesPair[j].second = j;
}
sort(eigenvaluesPair.begin(), eigenvaluesPair.end(), &comparator);
mypair tempPair = eigenvaluesPair[1];
for (unsigned int j=1; j<(N-1); j++) {
	eigenvaluesPair[j] = eigenvaluesPair[j+1];
}
eigenvaluesPair[N-1] = tempPair;
for (unsigned int j=0; j<N; j++) {
	eigenValues1[j] = eigenvaluesPair[j].first;
	eigenVectors1Complex.col(j) = eigenVectors1Copy.col(eigenvaluesPair[j].second);
}
//cout << "eigenValues1 = " << endl;
//cout << eigenValues1 << endl << endl;
for (unsigned int j=0; j<N; j++) {
	for (unsigned int k=0; k<N; k++) {
		eigenVectors1(j,k) = real(eigenVectors1Complex(j,k));
	}
}

//#pragma omp parallel for	
for (unsigned int j=0; j<N; j++) {
	for (unsigned int k=0; k<N; k++) {
		for (unsigned int l=0; l<N; l++) {
			double djdk = 4.0*pi*dr;
			if (j==0 || j==N) djdk/=sqrt(2.0);
			if (k==0 || k==N) djdk/=sqrt(2.0);
			omega1(j,k) += djdk*pow(eigenValues1(l),0.5)*eigenVectors1(j,l)*eigenVectors1(k,l);
			Eomega1(j,k) += djdk*eigenValues1(l)*eigenVectors1(j,l)*eigenVectors1(k,l);
		}
	}
}

/*---------------------------------------------------------------------------------------------
	2) analytic lattice expression based on ignoring b.c.s
---------------------------------------------------------------------------------------------*/

omega2 = Eigen::MatrixXd::Zero(N,N);
Eomega2 = Eigen::MatrixXd::Zero(N,N);
vec eigenValues2(N);
mat eigenVectors2(N,N);
double normalisation = sqrt(2.0/(N-1.0));
for (unsigned int l=0; l<N; l++) {
	eigenValues2(l) = 1.0+pow(2.0*sin(pi*l/(N-1.0)/2.0)/dr,2.0);
	for (unsigned int m=0; m<N; m++) eigenVectors2(l,m) = normalisation*sin(pi*l*m/(N-1.0));
}
eigenValues2(N-1) = 1.0;
//cout << "eigenValues2 = " << endl;
//cout << eigenValues2 << endl << endl;
for (unsigned int j=0; j<N; j++){
	for (unsigned int k=0; k<N; k++){
		for (unsigned int l=0; l<N; l++){
			double djdk = 4.0*pi*dr;
			omega2(j,k) += djdk*pow(eigenValues2(l),0.5)*eigenVectors2(j,l)*eigenVectors2(k,l);
			Eomega2(j,k) += djdk*eigenValues2(l)*eigenVectors2(j,l)*eigenVectors2(k,l);
		}
	}
}

/*---------------------------------------------------------------------------------------------
	checking norms
---------------------------------------------------------------------------------------------*/
vec norms1(N), norms2(N);
norms1 = eigenVectors1.rowwise().norm();
norms2 = eigenVectors2.rowwise().norm();
for (unsigned int j=0; j<N; j++) {
	if (absDiff(norms1(j),1.0)>1.0e-10) cout << "norms1(" << j << ") =  " << norms1(j) << endl;
	if (absDiff(norms2(j),1.0)>1.0e-10) cout << "norms2(" << j << ") =  " << norms2(j) << endl;
}
cout << endl;

/*---------------------------------------------------------------------------------------------
	calculating differences
---------------------------------------------------------------------------------------------*/

vec eigenvalueDiff(N);
for (unsigned int j=0; j<N; j++) {
	eigenvalueDiff(j) = absDiff(eigenValues1(j),eigenValues2(j));
	if (eigenvalueDiff(j)>1.0e-12) cout << "eigenvalueDiff(" << j << ") = " << eigenvalueDiff(j) << endl;
}
cout << endl;

vec eigenvectorDiff(N);
for (unsigned int j=0; j<N; j++) {
	eigenvectorDiff(j) = absDiff(eigenVectors1.row(j),eigenVectors2.row(j));
//	if (eigenvectorDiff(j)>1.0e-10) cout << "eigenvectorDiff(" << j << ") = " << eigenvectorDiff(j) << endl;
}
cout << endl;

vec eigenvectorSolution1(N);
for (unsigned int j=0; j<N; j++) {
	vec solutionDiff(N);
	vec eigenvectorTemp = eigenVectors1.col(j);
	solutionDiff = h*eigenvectorTemp;
	double approxEigenvalue = solutionDiff.dot(eigenvectorTemp);
	//solutionDiff -= eigenValues2(j)*eigenvectorTemp;
	solutionDiff -= approxEigenvalue*eigenvectorTemp;
	eigenvectorSolution1(j) = solutionDiff.norm();
	eigenvectorSolution1(j) /= 2.0*approxEigenvalue*N;
	if (eigenvectorSolution1(j)>1.0e-4) cout << "eigenvectorSolution1(" << j << ") = " << eigenvectorSolution1(j) << endl;
}
cout << endl;

vec eigenvectorSolution2(N);
for (unsigned int j=0; j<N; j++) {
	vec solutionDiff(N);
	vec eigenvectorTemp = eigenVectors2.col(j);
	solutionDiff = h*eigenvectorTemp;
	double approxEigenvalue = solutionDiff.dot(eigenvectorTemp);
	//solutionDiff -= eigenValues2(j)*eigenvectorTemp;
	solutionDiff -= approxEigenvalue*eigenvectorTemp;
	eigenvectorSolution2(j) = solutionDiff.norm();
	eigenvectorSolution2(j) /= 2.0*approxEigenvalue*N;
	if (eigenvectorSolution2(j)>1.0e-4) cout << "eigenvectorSolution2(" << j << ") = " << eigenvectorSolution2(j) << endl;
}
cout << endl;

mat omegaDiff = omega1-omega2;
Eigen::MatrixXd::Index maxRow, maxCol;
cout << "omegaDiff.maxCoeff() = " << omegaDiff.maxCoeff(&maxRow,&maxCol) << ", at (" << maxRow << "," << maxCol << ")" << endl;
cout << endl;

mat EomegaDiff = Eomega1-Eomega2;
cout << "EomegaDiff.maxCoeff() = " << EomegaDiff.maxCoeff(&maxRow,&maxCol) << ", at (" << maxRow << "," << maxCol << ")" << endl;
cout << endl;

// quick test to see column/row ordering
double quickTestRow, quickTestColumn;
vec quickTestRowVec(N), quickTestColumnVec(N), temp(N);
for (unsigned int j=0; j<N; j++) {
	quickTestRowVec(j) = eigenVectors1(2,j);
	quickTestColumnVec(j) = eigenVectors1(j,2);
}
temp = h*quickTestRowVec;
temp -= eigenValues1(2)*quickTestRowVec;
quickTestRow = temp.norm();
temp = h*quickTestColumnVec;
temp -= eigenValues1(2)*quickTestColumnVec;
quickTestColumn = temp.norm();
cout << "quickTestRow = " << quickTestRow << endl;
cout << "quickTestColumn = " << quickTestColumn << endl;
cout << endl;

//printing an eigenvalue of each
cout << "eigenValues1(0) = " << eigenValues1(0) << endl;
cout << "eigenValues2(0) = " << eigenValues2(0) << endl;
cout << endl;

return 0;
}
