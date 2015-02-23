/*
	simple program to read a file and output that which follows a given string
*/

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

void simpleReader (const string& fileIn, const string& fileOut, const int& nexts) {

}

int main() {

string fileIn = "results/13.02.15_pi_output.txt", fileOut = "temp", toSearch = "linErgA";
int nexts = 2;
/*
cout << "enter filename to search: ";
cin >> fileIn;*/
cout << "enter string to search for: ";
cin >> toSearch;/*
cout << "look how many inputs after string: ";
cin >> nexts;
cout << "enter filename to output to: ";
cin >> fileOut;
*/


ifstream is;
ofstream os;
is.open(fileIn.c_str());
os.open(fileOut.c_str());
if (!is.good()) cerr << fileIn << " not opened properly" << endl;
if (!os.good()) cerr << fileOut << " not opened properly" << endl;
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
