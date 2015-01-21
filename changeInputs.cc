// quick program to change the inputs file from bash

#include <iostream>
#include <fstream>
#include <vector>
#include "files.h"

using namespace std;

int main(int argc, char** argv) {
string filename, pName, pValue;

if (argc==3) {
	pName = argv[1];
	pValue = argv[2];
	filename = "inputs";
}
else if (argc>3) {
	for (int j=0; j<(int)(argc/2); j++) {
		string temp1 = argv[2*j+1];
		string temp2 = argv[2*j+2];
		if (temp1[0]=='-') temp1 = temp1.substr(1);
		if (temp1.compare("f")==0 || temp1.compare("file")==0) filename = temp2;
		else if (temp1.compare("v")==0 || temp1.compare("value")==0) pValue = temp2;
		else if (temp1.compare("n")==0 || temp1.compare("name")==0) pName = temp2;
		else {
			cerr << "input " << temp1 << " not understood" << endl;
			return 1;
		}
	}
}
else cerr << "inputs not understood" << endl;

if (pName[0]=='-') pName = pName.substr(1);
changeInputs("temp",pName,pValue,filename);
copyFile("temp",filename);

return 0;
}
