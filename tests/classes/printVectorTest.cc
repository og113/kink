// program to test whether one can print vectors as columns using <<

#include <fstream>
#include <iostream>
#include <iomanip>
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

ostream& operator<<(ostream& os, const vector< vector<double> >& vecVec) {
	uint size = vecVec[0].size();
	if (vecVec.size()>0) {
		for (uint j=1; j<vecVec.size(); j++) {
			if (vecVec[j].size()!=size) throw "print error";
		}
		for (uint k=0; k<size; k++) {
			for (uint j=0; j<vecVec.size(); j++) {
				os << setw(16) << vecVec[j][k];
			}
			os << endl;
		}
	}
	else {
		for (uint k=0; k<size; k++) {
			os << vecVec[0][k] << endl;
		}
	}
	return os;
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

vector< vector<double> > w(2);
w[0] = u;
w[1] = v;
cout << w;

return 0;
}
