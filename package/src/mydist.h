/*
my licensing information
*/

#ifndef MYDIST_H
#define MYDIST_H

double spDist(size_t m, double* Xa, double* Xb);

class sqDist
{
public:
	size_t nX, mX;
	double* X;
	// double minL;
	// constructor
	sqDist(SEXP sexpX);
	// sparse-matrix row pair-wise distances
	void row_spDist(size_t i, double* Li);
	// distance-matrix row pair-wise distances
	void row_d2Dist(size_t i, double* Li);
	// data-matrix row pair-wise distances
	void row_d1Dist(size_t i, double* Li);
};

#endif
