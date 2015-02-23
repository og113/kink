/*
	simple program to read a file and output that which follows a given string
	should probably write this in a more generic way so that it can be reused more easily
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

int simpleReader (const string& fileIn, const string& fileOut, const string& toSearch, const int& nexts) {
ifstream is;
ofstream os;
is.open(fileIn.c_str());
os.open(fileOut.c_str());
if (!is.good()) {
	cerr << fileIn << " not opened properly" << endl;
	return 1;
}
if (!os.good()) {
	cerr << fileOut << " not opened properly" << endl;
	return 1;
}
string line, temp;
while (getline(is,line)) {
		if(line[0] == '#') continue;
		if(line.empty()) continue;
		size_t searchPos = line.find(toSearch);
		if(searchPos==string::npos) continue;
		line = line.substr(searchPos+toSearch.size());
		istringstream ss(line);
		for (int j=0; j<nexts; j++) ss >> temp;
		os << temp << endl;
}

is.close();
os.close();

return 0;
}

int main() {

string date = "22.02.15";
string identifier = "cutoff";
string out_file = "results/" + date + "_" + identifier + "_assembled.dat";
int 	N 	= 130;
int 	Nb 	= 80;
double 	L 	= 5.0;
double 	Tb 	= 0.8;
string fileIn = "results/" + date + "_" + identifier + ".txt", action_file = "data/" + date + "_action_" + identifier + ".dat";
string tempOut1 = "data/" + date + "_" + identifier + "_tempOut1.txt", tempOut2 = "data/" + date + "_" + identifier + "_tempOut2.txt";
string tempOut3 = "data/" + date + "_" + identifier + "_tempOut3.txt";
bool is1Complex = true, is2Complex = true, is3Complex = false;
bool is3HalfAsRegular = true;

simpleReader(fileIn,tempOut1,"linErgAB",2);
simpleReader(fileIn,tempOut2,"linErg(0)",2);
simpleReader(fileIn,tempOut3,"CUTOFF",2);

ifstream is_1, is_2, is_3, is_action;
ofstream os;
is_1.open(tempOut1.c_str());
is_2.open(tempOut2.c_str());
is_3.open(tempOut3.c_str());
is_action.open(action_file.c_str());
os.open(out_file.c_str());
if (!is_1.good()) {
	cerr << tempOut1 << " not opened properly" << endl;
	return 1;
}
if (!is_2.good()) {
	cerr << tempOut2 << " not opened properly" << endl;
	return 1;
}
if (!is_3.good()) {
	cerr << tempOut3 << " not opened properly" << endl;
	return 1;
}
if (!is_action.good()) {
	cerr << action_file << " not opened properly" << endl;
	return 1;
}
if (!os.good()) {
	cerr << out_file << " not opened properly" << endl;
	return 1;
}
string line, action, temp1, temp2, temp3, dross;
os << left;
while (getline(is_action,line)) {
		istringstream ss(line);
		ss >> dross >> N >> dross >> L >> Tb >> dross >> dross >> dross >> dross >> action;
		is_1 >> temp1;
		is_2 >> temp2;
		if (is3HalfAsRegular) {
			is_3 >> temp3;
			is3HalfAsRegular = !is3HalfAsRegular;
		}
		else {
			is3HalfAsRegular = true;
		}
		if (is1Complex) temp1 = temp1.substr(1,temp1.find(",")-1);
		if (is2Complex) temp2 = temp2.substr(1,temp2.find(",")-1);
		if (is3Complex) temp3 = temp3.substr(1,temp3.find(",")-1);
		os << setw(16) << N << setw(16) << Nb << setw(16) << L << setw(16) << Tb;
		os << setw(16) << temp1 << setw(16)  << temp2 << setw(16) << temp3 << setw(16) << action << endl;
}

is_1.close();
is_2.close();
is_3.close();
is_action.close();
os.close();

return 0;
}
