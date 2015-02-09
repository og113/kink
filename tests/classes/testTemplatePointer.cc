/*
main for program to test interaction of templates and function pointers
*/

#include <iostream>
#include <string>
#include "templatePointer.h"
#include "fnptr.h"

using namespace std;

double (*fn) (double);

int main() {
int i;
double d;
string s;

fn = &function<double>;

i = function<int>(2);
d = fn(1.2);
s = function<string>("str");

cout << i << " " << d << " " << s << endl;

return 0;
}
