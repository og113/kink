#include <iostream>

using namespace std;

typedef unsigned int uint;

/*-------------------------------------------------------------------------------------------------------------------------
point class
-------------------------------------------------------------------------------------------------------------------------*/
class Point_2d {
public:
	Point_2d(): t(), x()	{}						// empty constructor
	Point_2d(const uint& t_in, const uint& x_in) {	// constructor
		t = t_in;
		x = x_in;
	}
	~Point_2d() {}									// empty destructor
	uint operator[](const uint& proj) const;		// projection operator
private:
	uint		t;									// coord values
	iint		x;									//
};

/*-------------------------------------------------------------------------------------------------------------------------
iterator over points
-------------------------------------------------------------------------------------------------------------------------*/
class Iterator_Point_2d {
public:
	Iterator_Point_2d(): Nt(), Nx()	{}				// empty constructor
	Iterator
private:
	uint		Nt;									// size of t dimension
	uint		Nx;									// size of x dimension
};
