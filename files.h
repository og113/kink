//functions related to files for pi.cc and main.cc
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <ctype.h>
#include <cctype>
#include <gsl/gsl_poly.h>
#include "gnuplot_i.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//simpler generic functions

//to convert number to string, usage is string str = NumberToString<number type>(x);
template <typename T>
string numberToString ( T Number )
	{
	stringstream ss;
	ss << Number;
	return ss.str();
	}

//to convert string to number, usage is (number type) x = StringToNumber<number type>(str);
template <typename T>
T stringToNumber ( const string &Text )//Text not by const reference so that the function can be used with a 
	{                               //character array as argument
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
	}

//getting the date and time
const string currentDateTime()
	{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%y%m%d%H%M%S", &tstruct);
    return buf;
	}
		
//copy a file
void copyFile(const string & inputFile, const string & outputFile)
	{
	ifstream  src(inputFile.c_str(), ios::binary);
	ofstream  dst(outputFile.c_str(), ios::binary);

	dst << src.rdbuf();
	}

//get last line of a file
string getLastLine(ifstream& inStream)
	{
    string xLine;
    while (inStream >> ws && getline(inStream, xLine)) // skip empty lines
    	{
        ;
		}		
    return xLine;
	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//file string manipulation

//get elements of string vector that have a given string in them
vector<string> findStrings(const vector <string> & fullVector, const string & search)
	{
	vector <string> subVector;
	for (unsigned int l=0;l<fullVector.size();l++)
		{
		if(fullVector[l].find(search)!=string::npos)
			{
			subVector.push_back(fullVector[l]);
			}
		}
	return subVector;
	}
	
//get vector of unsigned long long int from vector of strings
vector<unsigned long long int> getInts(const vector <string> & strVector)
	{
	vector <unsigned long long int> intVector;
	for (unsigned int l=0; l<strVector.size(); l++)
		{
		string temp = strVector[l];
		if (temp[7]=='1')
			{
			temp = temp.substr(7,12);
			intVector.push_back(stringToNumber<unsigned long long>(temp));
			}
		else
			{
			cout << "getInts error, filename not as expected: " << temp << endl;
			}
		}
	return intVector;
	}
	
//function to return final numbers in strings, after last "_", note it will also change argument
vector<int> getLastInts(vector <string> * strVector) //this function could be improved by splitting into two functions, one that finds the integers at the end of a string and the other that loops over the vector, ideally we wouldn't have to look for "_" but could just look for not digit, something along these lines would surely work
	{
	vector <int> intVector;
	for (unsigned int l=0; l<(*strVector).size(); l++)
		{
		size_t first_index;
		size_t last_index;
		string str = (*strVector)[l];
		if (str.find_last_of("_")!=string::npos)
			{
			first_index = str.find_last_of("_");
			}
		else
			{
			//cout << "getLastInt error, no underscore in file input so file removed from list:" << endl;
			//cout << (*strVector)[l] << endl;
			(*strVector)[l] = -1;
			continue;
			}
		str = str.substr(first_index + 1);
		if (isdigit(str[0])==0) //zero for false
			{
			//cout << "getLastInt error, character after _ is not digit so file removed from list:" << endl;
			//cout << (*strVector)[l] << endl;
			(*strVector)[l] = -1;
			continue;
			}
		if (str.find_last_of("0123456789")!=string::npos)
			{
			last_index = str.find_last_of("0123456789");
			}
		else
			{
			//cout << "getLastInt error, no numbers in file input so file removed from list:" << endl;
			//cout << (*strVector)[l] << endl;
			(*strVector)[l] = -1;
			continue;
			}
		str = str.substr(0,last_index+1);
		intVector.push_back(stringToNumber<unsigned long long>(str));		
		}
	return intVector;
	}
	
//function to reduce two vectors of strings to the size of another, only keeping elements with the timeNumber in common
vector<string> reduceTo(vector <string> toReduce, const vector <string> & toCompare)
	{
	vector <unsigned long long int> toReduceNumbers = getInts(toReduce);
	vector <unsigned long long int> toCompareNumbers = getInts(toCompare);
	if (toCompare.size()<toReduce.size() && toCompare.size()>0)
		{
		for (unsigned int j=0;j<toReduce.size();j++)
			{
			if(find(toCompare.begin(), toCompare.end(), toReduce[j]) == toCompare.end())
				{
				toReduce.erase(toReduce.begin()+j);
				}
			}
		}
	else if (toCompare.size()==0)
		{
		toReduce = toCompare;
		}
	return toReduce;
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//specific file manipulations
	
//read dataFiles into filenames and filenumbers
vector<string> readDataFiles(const unsigned long long int & minFileNo, const unsigned long long int & maxFileNo)
	{
	vector<string> fileNames;
	unsigned long long int fileNumber;
    ifstream file;
    file.open ("dataFiles");
    string fileName;
    string strNumber;
    fileName.clear();
		while ( !file.eof() )
			{
			file >> fileName;
			if (fileName.size()>19)
				{
				if (fileName[7]=='1' && fileName[fileName.size()-1]!='~')
					{
					strNumber = fileName.substr(7,12);
					fileNumber = stringToNumber<unsigned long long>(strNumber);
					if (fileNumber>=minFileNo && fileNumber<=maxFileNo)
						{
						fileNames.push_back(fileName);
						}
					}
				fileName.clear();
				}
    		}
    file.close();
    return fileNames;
	}
	
//copy inputs with a change
void changeInputs(const string & outputFile, const string & search, const string & replace)
	{ 
	ifstream fin;
	ofstream fout;
	fin.open("inputs");
	fout.open(outputFile.c_str());
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
	}
	

