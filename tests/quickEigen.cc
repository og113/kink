#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace std;

//factorial function
double factorial(const int& f_input)
	{
	double f_result = 1.0;
	for (int l=0; l<f_input; l++) f_result *= (double)(l+1.0);
	return f_result;
	}

//C^n_k
double Cnk(const int& n, const int& k) {
	if (n<k) cerr << "Ckn error: n=" << n << ", k=" << k << endl;
	return (double)factorial(n)/(double)factorial(k)/(double)factorial(n-k);
}

//comparators for complex<double>
bool operator<(const complex<double>& left, const complex<double>& right) {
	return real(left)<real(right);
}
bool operator>(const complex<double>& left, const complex<double>& right) {
	return real(left)>real(right);
}

// pair for sorting eigenvalues
typedef std::pair<complex<double>,size_t> mypair;
bool comparator ( const mypair& l, const mypair& r)
   { return l.first < r.first; }
   
//to convert number to string, usage is string str = NumberToString<number type>(x);
template <typename T>
string numberToString ( T Number ){
	stringstream ss;
	ss << Number;
	return ss.str();
}

//to convert string to number, usage is (number type) x = StringToNumber<number type>(str);
template <typename T>
T stringToNumber ( const string &Text ){
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
unsigned int N=5;
double mass2 = 1.0, dx = 0.5;
double a = pow(dx,-2.0), b = mass2 + 2.0*pow(dx,-2.0);
if (argc==2) N = atoi(argv[1]);
else if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("m")==0 || id.compare("mass2")==0) 	mass2 = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("a")==0) 						a = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("b")==0) 						b = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("N")==0)                        N = stringToNumber<unsigned int>(argv[2*j+2]);
		else if (id.compare("dx")==0) 						dx = stringToNumber<double>(argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}
a = pow(dx,-2.0), b = mass2 + 2.0*pow(dx,-2.0);
Eigen::MatrixXd A(N+2,N+2);
for (unsigned int j=0; j<(N+2); j++) {
	if (j==0) {
		A(j,j) = 1.0;
		A(j,j+1) = 0.0;
	}
	else if (j==(N+2-1)) {
		A(j,j) = 1.0;
		A(j,j-1) = 0.0;
	}
	else {
		A(j,j-1) = -a;
		A(j,j+1) = -a;
		A(j,j)   = b;
	}
}
cout << "Here is the " << N+2 << "x" << N+2 << " matrix, A:" << endl << A << endl << endl;
Eigen::EigenSolver<Eigen::MatrixXd> es(A);
Eigen::VectorXcd eigenvalues = es.eigenvalues();
Eigen::MatrixXcd eigenvectors(N+2,N+2), eigenvectorsCopy = es.eigenvectors();
vector<mypair> eigenvaluesPair(N+2);
for (unsigned int j=0; j<(N+2); j++) {
	eigenvaluesPair[j].first = eigenvalues[j];
	eigenvaluesPair[j].second = j;
}
sort(eigenvaluesPair.begin(), eigenvaluesPair.end(), &comparator);
for (unsigned int j=0; j<(N+2); j++) {
	eigenvalues[j] = eigenvaluesPair[j].first;
	eigenvectors.col(j) = eigenvectorsCopy.col(eigenvaluesPair[j].second);
}
cout << "The eigenvalues of A are:" << endl << eigenvalues << endl;
cout << "The matrix of eigenvectors, V, is:" << endl << eigenvectors << endl << endl;
cout << "The determinant is: " << A.determinant() << endl;
double determinant = 0.0;
if (N%2) {
	unsigned int k = (unsigned int)(N/2);
	for (unsigned int j=0; j<=k; j++) {
		determinant += Cnk(k+j+1,2*j+1)*pow(b,2.0*j+1.0)*pow(-pow(a,2.0),(double)(k-j));
	}
}
else {
	unsigned int k = N/2;
	for (unsigned int j=0; j<=k; j++) {
		determinant += Cnk(k+j,2*j)*pow(b,2.0*j)*pow(-pow(a,2.0),(double)(k-j));
	}
}
cout << "The analytic expression is: " << determinant << endl;
return 0;
}
