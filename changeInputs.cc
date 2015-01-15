// quick program to change the inputs file from bash

#include <iostream>
#include <fstream>
#include "files.h"

using namespace std;

int main(int argc, char** argv) {
string pName, pValue;
if (argc==3) pName = argv[1];
else {
	cerr << "must give 2 arguments to changeInputs: name and value" << endl;
	return 1;
}
pValue = argv[2];
changeInputs("temp",pName,pValue);
copyFile("temp","inputs");

return 0;
}
