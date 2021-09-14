/*
my licensing information
*/

#ifndef MYDIST_H
#define MYDIST_H

double spDist(unsigned int m, double* Xa, double* Xb);

class sqDist
{
public:

	unsigned int nX, mX;
	double* X;
	// constructor
	sqDist(SEXP sexpX);
	// destrucor
	~sqDist() {
		X = NULL;
	}

	// sparse-matrix row pair-wise distances
	void row_spDist(unsigned int i, double* Li);
	// distance-matrix row pair-wise distances
	void row_d2Dist(unsigned int i, double* Li);
	// data-matrix row pair-wise distances
	void row_d1Dist(unsigned int i, double* Li);
};

#endif
