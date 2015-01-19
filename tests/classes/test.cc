// main program to do class tests

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "bag.h"

/*-------------------------------------------------------------------------------------------------------------------------
main
-------------------------------------------------------------------------------------------------------------------------*/

int main(){
Parameter<double> p("csi",1.0);
Parameter<double> q(p);
Parameter<double> r;
r = q;
cout << p << q << r;
Parameter<double> s("chi",2.0), t("pi",3.0), u("csi",0.0);
ParameterBag<double> b, c;
b.set(s);
b.set(t);
cout << b.size() << endl;
cout << b;
b.set(p);
b.set(u);
cout << b.size() << endl;
cout << b;
try{
cout << b("csi") << endl;
} catch (BagError::NotFound e) {
cerr << e;
}
try{cout << b("fi") << endl;
} catch (BagError::NotFound e) {
cerr << e;
}
b.reset();
b.set(r);
c = b;
cout << c;
c.reset();

ifstream file;
file.open("inputs");
if (file.good()) {
	file >> c;
	file.close();
}
else cerr << "inputs doesn't exist" << endl;
cout << c;

ofstream out;
out.open("outputs");
out << c;
out.close();


return 0;
}

