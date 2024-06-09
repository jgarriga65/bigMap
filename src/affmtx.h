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
	unsigned int z, nX, mX, nnSize;
	int* zIdx;
	double* X;
	double* B;
	double zP;
	// debug index
	int debugRow;

	// constructor
	affMtx(SEXP sexpX, SEXP sexpB, int* zIdx, unsigned int z, unsigned int nnSize);
	// destructor
	~affMtx() {
		zIdx = NULL;
		B = NULL;
		X = NULL;
	}

	// transform input similarities into unnormalized probabilities from INPUT-DATA
	void X2P(double* P, unsigned int* W);
	// transform input similarities into unnormalized probabilities from DISTANCES Triangular-Matrix
	void D2P(double* P, unsigned int* W);
	// transform input similarities into unnormalized probabilities from SPARSE-DATA
	void S2P(double* P, unsigned int* W);
	// affinities
	void affinities(double* K, double* P, unsigned int* W);
	// versions for BHt-SNE
	void bh_X2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W);
	void bh_D2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W);
	void bh_S2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W);

};

#endif
