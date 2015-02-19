// quick program to compare two vectors loaded from files

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

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

int main(int argc, char** argv){

string fileA = "data/0pip_0.dat", fileB = "data/1pip_0.dat";
int lengthA, lengthB;
double A, B, C, D, result = 0.0, norm = 0.0;
int colA = 4, colB = 4;
bool complex = true;

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
			else if (temp1.compare("colA")==0) colA = stringToNumber<int>(temp2);
			else if (temp1.compare("colB")==0) colB = stringToNumber<int>(temp2);
			else if (temp1.compare("c")==0) complex = (bool)stringToNumber<int>(temp2);
			else {
				cerr << "input " << temp1 << " not understood" << endl;
				return 1;
			}
		}
	}
}

lengthA = countLines(fileA);
lengthB = countLines(fileB);

if (lengthA!=lengthB) {
	cerr << "files not the same length: lengthA = " << lengthA << ", lengthB = " << lengthB << endl;
	return 1;
}

ifstream isA, isB;
string lineA, lineB;
isA.open(fileA.c_str());
isB.open(fileB.c_str());
while(!isA.eof()) {
		string temp;
		for (int j=0; j<colA; j++) isA >> temp;
		isA >> A;
		if (complex) isA >> C;
		for (int j=0; j<colB; j++) isB >> temp;
		isB >> B;
		if (complex) isB >> D;
		result += pow(A-B,2.0);
		norm += pow(A+B,2.0);
		if (complex) {
			result += pow(C-D,2.0);
			norm += pow(C+D,2.0);
		}
	}
result = sqrt(result);
norm = sqrt(norm);
result /= norm;
result *= 2.0;

cout << "2*norm(A-B)/norm(A+B) = " << result << endl;

}
