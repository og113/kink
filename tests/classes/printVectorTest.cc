// program to test whether one can print vectors as columns using <<

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

typedef unsigned int uint;

template< typename T, size_t N >
vector<T> makeVector( const T (&data)[N] )
{
    return vector<T>(data, data+N);
}

class myVec{
public:
	myVec(vector<double> vec) : the_vec(vec) {}
	~myVec() {}
	vector<double> the_vec;
};

ostream& operator<<(ostream& os, const myVec& vec) {
	for (uint j=0; j<vec.the_vec.size(); j++) {
		os << vec.the_vec[j] << endl;
	}
}

int main(){

const double tmpu[] = { 4.0, 5.0, 6.0 };
const double tmpv[] = { 1.0, 2.0, 3.0 };
vector<double> v = makeVector(tmpv);
vector<double> u = makeVector(tmpu);
myVec myv(v);
myVec myu(u);
cout << myv;
ofstream f;
f.open("outputs");
f << myv;
f.close();

return 0;
}
