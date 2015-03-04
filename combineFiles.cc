/*
	simple program to read two files and combine them
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

int main() {
string fileIn1 = "data/150223114037mainb_kE_0_0_7.dat", fileIn2 = "data/150223114037maina_kE_0_0_7.dat";
string fileOut = "results/23.02.15_cutoff_ab_2.dat";

ifstream is1, is2;
ofstream os;
is1.open(fileIn1.c_str());
is2.open(fileIn2.c_str());
os.open(fileOut.c_str());
if (!is1.good()) {
	cerr << fileIn1 << " not opened properly" << endl;
	return 1;
}
if (!is2.good()) {
	cerr << fileIn2 << " not opened properly" << endl;
	return 1;
}

if (!os.good()) {
	cerr << fileOut << " not opened properly" << endl;
	return 1;
}

string line1, line2;
while(!is1.eof() && !is2.eof()) {
	getline(is1,line1);
	getline(is2,line2);
	os << line1 << line2 << endl;
}

is1.close();
is2.close();
os.close();

return 0;
}
