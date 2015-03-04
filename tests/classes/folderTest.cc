/*
	program to test the classes and functions defined in folder.h and folder.cc
*/

#include <iostream>
#include <vector>
#include <utility>
#include "folder.h"

using namespace std;

int main() {
cout << "folderTest.cc:" << endl;
Filename f;
f.Directory = "data";
f.Timenumber = "";
f.ID = "piE";
vector<StringPair> exs(1);
exs[0] = StringPair("N","100");
f.Extras = exs;
f.Suffix = ".dat";
cout << f << endl;
Filename g(f);
g.ID = "minusDSE";
cout << g << endl;
Filename h(f());
h.Timenumber = "0";
cout << h << endl;
return 0;
}
