#include <iostream>

using namespace std;

uint Point_2d::operator[](const uint& proj) {
	if (proj>1) {
		cerr << "Point_2d error: projected dimension > 1" << endl;
		return 1;
	}
	(proj==0)? return t: return x;
}
