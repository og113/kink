/*
	simple program to read a file and output that which follows a given string
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

string out_file = "results/13.02.15_pi_assembled.dat";
int 	N 	= 130;
int 	Nb 	= 80;
double 	L 	= 5.0;
double 	Tb 	= 0.8;
string erg_file = "data/13.02.15_linErgA.txt", num_file = "data/13.02.15_linNumA.txt", action_file = "data/action.dat";

simpleReader("results/13.02.15_pi_output.txt",erg_file,"linErgA",2);
simpleReader("results/13.02.15_pi_output.txt",num_file,"linNumA",2);


ifstream is_erg, is_num, is_action;
ofstream os;
is_erg.open(erg_file.c_str());
is_num.open(num_file.c_str());
is_action.open(action_file.c_str());
os.open(out_file.c_str());
if (!is_erg.good()) {
	cerr << erg_file << " not opened properly" << endl;
	return 1;
}
if (!is_num.good()) {
	cerr << num_file << " not opened properly" << endl;
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
string line, action, erg, num, temp;
os << left;
while (getline(is_action,line)) {
		istringstream ss(line);
		ss >> temp >> N >> temp >> L >> Tb >> temp >> temp >> action;
		is_erg >> erg;
		is_num >> num;
		os << setw(16) << N << setw(16) << Nb << setw(16) << L << setw(16) << Tb;
		os << setw(16) << erg << setw(16)  << num << setw(16) << action << endl;
		Tb += 0.005;
}

is_erg.close();
is_num.close();
is_action.close();
os.close();

return 0;
}
