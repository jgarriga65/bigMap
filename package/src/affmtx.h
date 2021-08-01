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

// // sort neighours by affinity
// static inline std::vector<int> row_affinity_sort(size_t z, size_t nSS, double* Ki)
// {
// 	int r = 0;
// 	std::vector<int> rank(z);
// 	std::iota(rank.begin(), rank.end(), r ++);
// 	std::nth_element(rank.begin(), rank.begin() +nSS, rank.end(),
// 			[&](int a, int b) {return Ki[a] > Ki[b];}
// 	);
// 	return rank;
// }
// 
// static inline void row_affinity(unsigned int z, unsigned int nnSize, unsigned int zi, double* Ki, double* P, unsigned int* W, unsigned int &w, double &zP)
// {
// 	//sort
// 	std::vector<int> rank = row_affinity_sort(z, nnSize, Ki);
// 	// select and normalize
// 	for (size_t k = 0; ((k < nnSize) && (Ki[rank[k]] > .0)); k++) {
// 		unsigned int zj = rank[k];
// 		unsigned int ij = ijIdx(z, zi, rank[k]);
// 		if (P[ij] == .0) {
// 			W[w] = (zi < zj) ? zi *z +zj : zj *z +zi;
// 			w ++;
// 		}
// 		P[ij] += Ki[rank[k]];
// 		zP += P[ij];
// 	}
// }

class affMtx
{
public:

	// t-SNE input/ouput matrix size
	unsigned int z, w, nX, mX, nnSize;
	int* zIdx;
	double* X;
	double* B;
	double zP;

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
	// versions for BHt-SNE
	void bh_X2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W);
	void bh_D2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W);
	void bh_S2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W);

};

#endif
