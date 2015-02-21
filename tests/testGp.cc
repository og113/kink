//parameters and functions for pi.cc
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "gnuplot_i.hpp"

using namespace std;

int main()
{
string readFile = "./data/testErg.dat";
string commandOpenStr = "gnuplot -persistent";
const char * commandOpen = commandOpenStr.c_str();
FILE * gnuplotPipe = popen (commandOpen,"w");
string command1Str = "plot \"" + readFile + "\" with lines";
string command2Str = "pause -1";
const char * command1 = command1Str.c_str();
const char * command2 = command2Str.c_str();
fprintf(gnuplotPipe, "%s \n", command1);
fprintf(gnuplotPipe, "%s \n", command2);
pclose(gnuplotPipe);
return 0;
}
