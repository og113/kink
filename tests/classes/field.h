/*
classes and functions for the implementation of a 2d complex field, also containing zero mode lagrange multipliers, using the eigen vector class
*/

#ifndef __FIELD_H_INCLUDED__
#define __FIELD_H_INCLUDED__

#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include "error.h"

using namespace std;

typedef unsigned int uint;
typedef complex<double> comp;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;

/*-------------------------------------------------------------------------------------------------------------------------
simple coord functions
-------------------------------------------------------------------------------------------------------------------------*/

uint invCoord(const uint& t, const uint& x, const uint& Nt) {
	return t + x*Nt;
}

/*-------------------------------------------------------------------------------------------------------------------------
field error
-------------------------------------------------------------------------------------------------------------------------*/

class FieldError {
public:
	class OutOfRange: public SimpleError{
	public:
		OutOfRange(const string& rangeName, const uint& index, const uint& range)
			: the_rangeName(rangeName), the_index(index), the_range(range) {}	//constructor
		virtual string		message() const;									// message to be passed for printing
	private:
		string	the_rangeName;
		uint	the_index;
		uint	the_range;
	};
};


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
	uint			size() const;							// getting size of field vector
	uint			Nt()		{return the_Nt;}	const	// getting Nt
	uint			Nx()		{return the_Nx;}	const	// getting Nx
	uint			extras()	{return the_extras;}	const// getting extras
	vec&			vector()	{return &fieldVector;}		// getting address of fieldVector
private:
	uint			the_Nt;							// number of points in t direction
	uint			the_Nx;							// number of points in x direction
	uint			the_extras;						// number of extra field values, e.g. lagrange multipliers
	vec				fieldVector;					// the real and imagninary components of the field
};

ostream& operator<<(ostream& os, const Field& f);			// printing field
istream& operator<<(istream& is, Field& f);					// getting field values from stream

#endif // __FIELD_H_INCLUDED__
