// quick program to compare two sparse matrices loaded from files

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Sparse>

using namespace std;

typedef unsigned int uint;

//count non-empty lines of a file
unsigned int countLines(const string & file_to_count) {
	ifstream fin;
	fin.open(file_to_count.c_str());
	if (!fin.good()) {
		cerr << "countLines error: file " << file_to_count << " not opened properly" << endl;
		return 0;
	}
	string line;
	unsigned int counter = 0;
	while(!fin.eof()) {
		getline(fin,line);
		if(line.empty()) continue;
		counter++;
	}		
	fin.close();
    return counter;
}

//to convert string to number, usage is (number type) x = StringToNumber<number type>(str);
template <typename T>
T stringToNumber ( const string &Text )//Text not by const reference so that the function can be used with a 
	{                               //character array as argument
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
	}



/*------------------------------------------------------------------------------------------------------------------------
main begins
------------------------------------------------------------------------------------------------------------------------*/
int main(int argc, char** argv){

string fileA = "data/0DDSp_0.dat", fileB = "data/1DDSp_0.dat";
unsigned int lengthA = 0, lengthB = 0, i, j;
double result, normA, normB, v;
double closeness = 1.0e-16;

if (argc>1) {
	if (argc==3) {
	fileA = argv[1];
	fileB = argv[2];
	}
	else {
		for (int j=0; j<(int)(argc/2); j++) {
			string temp1 = argv[2*j+1];
			string temp2 = argv[2*j+2];
			if (temp1[0]=='-') temp1 = temp1.substr(1);
			if (temp1.compare("fileA")==0) fileA = temp2;
			else if (temp1.compare("fileB")==0) fileB = temp2;
			else {
				cerr << "input " << temp1 << " not understood" << endl;
				return 1;
			}
		}
	}
}

ifstream isA;
isA.open(fileA.c_str());
if (!isA.good()) {
	cerr << fileA << " not opened properly" << endl;
	return 1;
}
string line;
while(getline(isA,line)) {
		if(line.empty()) continue;
		istringstream ss(line);
		ss >> i >> j;
		if (i>lengthA) lengthA = i;
		if (j>lengthA) lengthA = j;
	}
isA.close();

ifstream isB;
isA.open(fileB.c_str());
if (!isB.good()) {
	cerr << fileB << " not opened properly" << endl;
	return 1;
}
while(getline(isB,line)) {
		if(line.empty()) continue;
		istringstream ss(line);
		ss >> i >> j;
		if (i>lengthB) lengthB = i;
		if (j>lengthB) lengthB = j;
	}
isB.close();

if (lengthA!=lengthB) {
	cerr << "lengthA(" << lengthA << ") != lengthB(" << lengthB << ")" << endl;
	return 1;
}

Eigen::VectorXi A_to_reserve = Eigen::VectorXi::Zero(lengthA);
Eigen::VectorXi B_to_reserve = Eigen::VectorXi::Zero(lengthA);
Eigen::VectorXi diff_to_reserve = Eigen::VectorXi::Zero(lengthA);

isA.open(fileA.c_str());
while(getline(isA,line)) {
		if(line.empty()) continue;
		istringstream ss(line);
		ss >> i;
		A_to_reserve(i)++;
	}
isA.close();

isB.open(fileB.c_str());
while(getline(isB,line)) {
		if(line.empty()) continue;
		istringstream ss(line);
		ss >> i;
		B_to_reserve(i)++;
	}
isB.close();

for (unsigned int j=0; j<lengthA; j++) diff_to_reserve(j) = A_to_reserve(j) + B_to_reserve(j);

Eigen::SparseMatrix<double> matA(lengthA,lengthA), matB(lengthA,lengthA), diff(lengthA,lengthA);

matA.reserve(A_to_reserve);
matB.reserve(B_to_reserve);
diff.reserve(diff_to_reserve);

string lineA, lineB;
isA.open(fileA.c_str());
isB.open(fileB.c_str());
while(!isA.eof()) {
		isA >> i >> j >> v;
		matA.insert(i,j) = v;
	}
while(!isB.eof()) {
		isB >> i >> j >> v;
		matB.insert(i,j) = v;
	}
isA.close();
isB.close();

matA.makeCompressed();
matA.prune(closeness);
matB.makeCompressed();
matB.prune(closeness);

diff = matA-matB;
diff.makeCompressed();
diff.prune(closeness);

normA = matA.squaredNorm();
normB = matB.squaredNorm();
result = diff.squaredNorm();
result /= sqrt(normA+normB);
result *= 2.0;

cout << "matA.nonZeros() = " << matA.nonZeros() << endl;
cout << "matB.nonZeros() = " << matB.nonZeros() << endl;
cout << "diff.nonZeros() = " << diff.nonZeros() << endl;
cout << "2*norm(A-B)/sqrt(norm(A)^2+norm(B)^2) = " << result << endl;

}
