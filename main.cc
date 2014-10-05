#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "pf.h"

int main()
{
ifstream fmainin;
fmainin.open("mainInputs", ios::in);
string line;
unsigned long long int minFile;
unsigned long long int maxFile;
while(getline(fmainin,line))
	{
	if(line[0] == '#')
		{
		continue;
		}
	if (!line.empty())
		{
		istringstream ss(line);
		ss >> minFile >> maxFile;
		}
	}
fmainin.close();


system("dir ./data/* > dataFiles");
vector<string> filenames, piFiles, inputsFiles;
filenames = readDataFiles(minFile,maxFile);
piFiles = findStrings(filenames,"tpip");
inputsFiles = findStrings(filenames,"inputs");

return 0;
}

