/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

#ifndef AFFMTX_H
#define AFFMTX_H

static inline unsigned int ijIdx(unsigned int z, unsigned int i, unsigned int j)
{
	return (i < j) ? (z *i) -(i +1) *(i +2) /2 +j : (z *j) -(j +1) *(j +2) /2 +i ;
}

class affMtx
{
public:

	// t-SNE input/ouput matrix size
	unsigned int z, w, nX, mX, nnSize;
	double latEx;
	int* zIdx;
	double* X;
	double* B;
	double zP;

	// constructor
	affMtx(SEXP sexpX, SEXP sexpB, int* zIdx, unsigned int z, unsigned int nnSize, double latEx);
	// destructor
	~affMtx() {
		zIdx = NULL;
		B = NULL;
		X = NULL;
	}

	// exaggeration factor
	void exgg(double* E);
	// transform input similarities into unnormalized probabilities from INPUT-DATA
	void X2P(double* P, int* W);
	// transform input similarities into unnormalized probabilities from DISTANCES Triangular-Matrix
	void D2P(double* P, int* W);
	// transform input similarities into unnormalized probabilities from SPARSE-DATA
	void S2P(double* P, int* W);
	// versions for enbedding final compression
	void efr_X2P(unsigned int z_ini, unsigned int z_end, double* P, int* W);
	void efr_D2P(unsigned int z_ini, unsigned int z_end, double* P, int* W);
	void efr_S2P(unsigned int z_ini, unsigned int z_end, double* P, int* W);

};

#endif
