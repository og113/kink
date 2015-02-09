/*
	quick practice for std::map
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <map>

using namespace std;

template <class Key, class Value>
ostream& operator<<(ostream& os, pair<Key,Value> p) {
	cout << left << setw(20) << p.first << setw(20) << p.second;
	return os;
} 

int main() {
map<string,double> m;
m.insert (pair<string,double>("a",0.1));

for (map<string,double>::iterator it=m.begin(); it!=m.end(); ++it)
	cout << *it << endl;

string a = "123456123456", b = "123456123457", c = "123456123455";

if (a>b) cout << "a>b" << endl;
else if (b>a) cout << "b>a" << endl;
else cout << "a==b" << endl;
if (a>c) cout << "a>c" << endl;
else if (c>a) cout << "c>a" << endl;
else cout << "a==c" << endl;

return 0;
}
