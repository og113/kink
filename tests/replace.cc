//test function to read inputs file and replace something
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main()
{
string search = "theta";
string replace = "0.5";

ifstream fin;
ofstream fout;
fin.open("inputs");
fout.open("./data/testInputs");
string line;
size_t pos;
while(!fin.eof())
	{
	getline(fin,line);
	if(line[0] == '#' && line[1] == '#')
		{
		fout << line << endl;
		continue;
		}
	if (line[0] == '#')
		{
		pos = line.find(search);
		fout << line << endl;
		getline(fin,line);
		if (pos != string::npos)
			{
			line.replace(pos, search.length(), replace);
			}
		}
	fout << line << endl;
	}
fin.close();
fout.close();

return 0;
}
