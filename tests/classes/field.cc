// implementation of the field class and related classes and functions

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <Eigen/Dense>
#include "simple.h"
#include "field.h"

complex<double> ii(0.0,1.0);

/*-------------------------------------------------------------------------------------------------------------------------
field error
-------------------------------------------------------------------------------------------------------------------------*/

FieldError::OutOfRange::message() {
	string s = "Field index out of " + the_rangeName + " range: index = " + NumberToString<uint>(index) \
				+ ", range = " + NumberToString<uint>(range);
	return s;
}

/*-------------------------------------------------------------------------------------------------------------------------
field
-------------------------------------------------------------------------------------------------------------------------*/

class Field {
public:
	Field(const uint& Nt_in, const uint& Nx_in, const uint& extras_in)
		:	the_Nt(Nt_in), the_Nx(Nx_in), the_extras(extras_in), fieldVector(2*Nt_in*Nx_in+extras) {}	// constructor
	~Field() {~fieldVector()}								// destructor
	operator comp	()(const uint& t, const uint& x) const;	// access complex component
	operator comp&	()(const uint& t, const uint& x);		//
	double 			r(const uint& t, const uint& x) const;	// access real part of component
	double&			r(const uint& t, const uint& x);		//
	double 			i(const uint& t, const uint& x) const;	// access imaginary part of component
	double&			i(const uint& t, const uint& x);		//
	double			extra(const uint& m) const;				// getting the extra components
	double&			extra(const uint& m);					//
	void			zero();									// sets all components to zero
	uint			size()		{return fieldVector.size()} const	// getting size of field vector
	uint			Nt()		{return the_Nt;}	const	// getting Nt
	uint			Nx()		{return the_Nx;}	const	// getting Nx
	uint			extras()	{return the_extras;}	const// getting extras
	vec&			vector()	{return &fieldVector;}		// getting address of fieldVector
	ostream& operator<<(ostream& os, const Field& f)	const;// printing field
	istream& operator<<(istream& is, Field& f);				// getting field values from stream
private:
	uint			the_Nt;							// number of points in t direction
	uint			the_Nx;							// number of points in x direction
	uint			the_extras;						// number of extra field values, e.g. lagrange multipliers
	vec				fieldVector;					// the real and imagninary components of the field
};

comp Field::operator()(const uint& t, const uint& x) const {
	if (t>the_Nt) throw FieldError::OutOfRange("t",t,the_Nt);
	if (x>the_Nx) throw FieldError::OutOfRange("x",x,the_Nx);
	uint pos = invCoord(t,x,the_Nt);
	return fieldVector(2*pos) + ii*fieldVector(2*pos+1);
}

double Field::r(const uint& t, const uint& x) const {
	if (t>the_Nt) throw FieldError::OutOfRange("t",t,the_Nt);
	if (x>the_Nx) throw FieldError::OutOfRange("x",x,the_Nx);
	uint pos = invCoord(t,x,the_Nt);
	return fieldVector(2*pos);
}

double Field::i(const uint& t, const uint& x) const {
	if (t>the_Nt) throw FieldError::OutOfRange("t",t,the_Nt);
	if (x>the_Nx) throw FieldError::OutOfRange("x",x,the_Nx);
	uint pos = invCoord(t,x,the_Nt);
	return fieldVector(2*pos+1);
}

double Field::extra(const uint& m) const {
	if (m>the_extras) throw FieldError::OutOfRange("extras",m,the_extras);
	return fieldVector(2*the_Nt*the_Nx+m);
}

void Field::zero() {
	fieldVector = Eigen::VectorXd::Zero(2*the_Nt*the_Nx+the_extras);
}

uint Field::size() const {
	return fieldVector.size();
}

ostream& operator<<(ostream& os, const Field& f) {

}

istream& operator<<(istream& is, Field& f) {

}
