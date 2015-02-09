/*
	quick test of map_extras
*/

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "map_extras.h"

struct A {
	double a;
};

struct B {
	double b;
};

struct C: A, B {
	C(): A(), B() {}
	C(A& Ast,B& Bst): A(Ast), B(Bst) {}
};

int main() {

map<string,double> m, n, r, s;
pair<string,double> p("a",0.1), q("b",0.2);

m.insert(p);
m.insert(q);
n.insert(q);

string file = "out.txt";

ofstream fout;
fout.open(file.c_str());
if (!fout.good()) cerr << "out.txt not opened properly" << endl;
else fout << m << n;
fout.close();

ifstream fin;
fin.open("out.txt");
if (!fin.good()) cerr << "out.txt not opened properly" << endl;
fin >> r;
fin >> s;
fin.close();

cout << r << s;

struct A Astr;
Astr.a=0.1;
struct B Bstr;
Bstr.b = 1.0;
struct C Cstr(Astr,Bstr);
cout << Cstr.a << " " << Cstr.b << endl;

struct C Cstr2(Cstr);
cout << Cstr2.a << " " << Cstr2.b << endl;

struct C Cstr3;
Cstr3 = Cstr2;
cout << Cstr3.a << " " << Cstr3.b << endl;
return 0;
}
