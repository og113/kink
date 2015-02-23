/*
	program to quickly test random crap
*/

#include <iostream>

using namespace std;

int main() {

unsigned int Na = 320;
unsigned int linearInt = (unsigned int)(Na/10);
cout << linearInt << endl;
double E_exact = 0.0;
for (unsigned int j=1; j<(linearInt+1); j++) E_exact += 0.1;
cout << E_exact << endl;
E_exact /= (double)linearInt;
cout << E_exact << endl;

return 0;
}
