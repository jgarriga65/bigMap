/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <bigmemory/MatrixAccessor.hpp>
#include "bigmemory/BigMatrix.h"
// [[Rcpp::depends(BH, bigmemory)]]

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <math.h>
#include <float.h>
#include "sqdist.h"
#include "affmtx.h"

using namespace Rcpp;


// affMtx constructor
affMtx::affMtx(SEXP sexpX, SEXP sexpB, int* zIdx, unsigned int z, unsigned int nnSize): zIdx(zIdx), z(z), nnSize(nnSize)
{

	// Att!! X and B transposed
	index_type offset;

	BigMatrix *bigX = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpX));
	offset = bigX ->nrow() *bigX ->col_offset();
	this ->X = reinterpret_cast <double*> (bigX ->matrix()) +offset;
	this ->nX = bigX ->ncol();
	this ->mX = bigX ->nrow();

	BigMatrix *bigB = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpB));
	offset = bigB ->nrow() *bigB ->col_offset();
	this ->B = reinterpret_cast <double*> (bigB ->matrix()) +offset;

	// thread affinities normalization factor
	this ->zP = .0;

	// clean
	bigX = NULL;
	bigB = NULL;
}

void affMtx::affinities(double* K, double* P, unsigned int* W)
{
	// row-affinities normalization factor
	std::vector<double> C(z, .0);
	// row-affinities
	for (unsigned int zi = 0; zi < z; zi++) {
		// sort neighbours by affinity
		// (bi-direccional 2)
		// std::vector<double> Aff(z);
		// for (unsigned int zj = 0; zj < z; zj++) Aff[zj] = K[zi *z +zj] + K[zj *z +zi];
		int r = 0;
		std::vector<int> rank(z);
		std::iota(rank.begin(), rank.end(), r++);
		std::nth_element(rank.begin(), rank.begin() +nnSize, rank.end(),
		//  (uni-directional)
		// 	[&](int a, int b) {return K[zi *z +a] > K[zi *z +b];}
		//  (bi-directional)
			[&](int a, int b) {return K[zi *z +a] +K[a *z +zi] > K[zi *z +b] +K[b *z +zi];}
		//	(bi-directional 2)
		// 	[&Aff](int a, int b) {return Aff[a] > Aff[b];}
		);
		// select most affine neighbours and compute normalization factor
		for (unsigned int ni = 0, ij = zi *nnSize; ni < nnSize; ni++, ij++) {
			W[ij] = rank[ni];
			C[zi] += K[zi *z +rank[ni]];
		}
	}
	double joinNormalization = .5 /z; 
	for (unsigned int zi = 0; zi < z; zi++) {
		for (unsigned int ni = 0, ij = zi *nnSize; ni < nnSize; ni++, ij++) {
			unsigned int zj = W[ij];
			P[ij] = joinNormalization *(K[zi *z +zj] /C[zi] + K[zj *z +zi] /C[zj]);
			zP += P[ij];
		}
	}
}

// +++ exact NN
void affMtx::X2P(double* P, unsigned int* W)
{
	std::vector<double> K(z *z, .0);
	for (unsigned int zi = 0; zi < z; zi++) {
		double Bi = B[zIdx[zi] *3 +0];
		double Ci = B[zIdx[zi] *3 +1];
		double Xi[mX];
		for(unsigned int v = 0, mi = zIdx[zi] *mX; v < mX; v++) Xi[v] = X[mi +v];
		for (unsigned int zj = zi +1; zj < z; zj++) {
			double Lij = .0;
			for(unsigned int v = 0, mj = zIdx[zj] *mX; v < mX; v++) Lij += pow(Xi[v] -X[mj +v], 2);
			K[zi *z +zj] = std::exp(-Bi *Lij) /Ci;
			K[zj *z +zi] = std::exp(-B[zIdx[zj] *3 +0] *Lij) /B[zIdx[zj] *3 +1];
		}
	}
	affMtx::affinities(K.data(), P, W);
}

// transform input euclidean-distances into probabilities
// FROM SPARSE-MATRIX DATA
void affMtx::S2P(double* P, unsigned int* W)
{
	std::vector<double> K(z *z, .0);
	for (unsigned int zi = 0; zi < z; zi++) {
		double Bi = B[zIdx[zi] *3];
		double Ci = B[zIdx[zi] *3 +1];
		double Xi[mX];
		for(unsigned int v = 0, mi = zIdx[zi] *mX; v < mX; v++) Xi[v] = X[mi +v];
		for (unsigned int zj = zi +1; zj < z; zj++) {
			double Xj[mX];
			for(unsigned int v = 0, mj = zIdx[zj] *mX; v < mX; v++) Xj[v] = X[mj +v];
			double Lij = spDist(mX, Xi, Xj);
			K[zi *z +zj] = std::exp(-Bi *Lij) /Ci;
			K[zj *z +zi] = std::exp(-B[zIdx[zj] *3 +0] *Lij) /B[zIdx[zj] *3 +1];
		}
	}
	affMtx::affinities(K.data(), P, W);
}

// transform input similarities into probabilities
// FROM FULL-DISTANCE-MATRIX
void affMtx::D2P(double* P, unsigned int* W)
{
	std::vector<double> K(z *z, .0);
	for (unsigned int zi = 0; zi < z; zi++) {
		double Bi = B[zIdx[zi] *3 +0];
		double Ci = B[zIdx[zi] *3 +1];
		for (unsigned int zj = zi +1; zj < z; zj++) {
			double Lij = std::pow(X[zIdx[zi] *mX +zIdx[zj]], 2);
			K[zi *z +zj] = std::exp(-Bi *Lij) /Ci;
			K[zj *z +zi] = std::exp(-B[zIdx[zj] *3 +0] *Lij) /B[zIdx[zj] *3 +1];
		}
	}
	affMtx::affinities(K.data(), P, W);
}

// +++ mt-SNE affinities +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// transform input similarities into probabilities
// FROM INPUT-DATA
void affMtx::mt_X2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W)
{
	for (unsigned int i = z_ini, zi = 0; i < z_end; i++, zi++) {
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[i *mX +v];
		double Bi = B[i *3 +0];
		double Ci = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int j = 0, ni = 0; ((j < nX) && (ni < nnSize)); j++) {
			if (j != i) {
				double Lij = .0;
				for(unsigned int v = 0; v < mX; v++) Lij += std::pow(Xi[v] -X[j *mX +v], 2);
				if (Lij <= Li) {
					unsigned int ij = zi *nnSize +ni;
					P[ij] += std::exp(-Bi *Lij) /Ci;
					P[ij] += std::exp(-B[j *3 +0] *Lij) /B[j *3 +1];
					zP += P[ij];
					W[ij] = j;
					ni++;
				}
			}
		}
	}
}


// transform input similarities into probabilities
// FROM FULL-DISTANCE-MATRIX
void affMtx::mt_D2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W)
{
	for (unsigned int i = z_ini, zi = 0; i < z_end; i++, zi++) {
		double Bi = B[i *3 +0];
		double Ci = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int j = 0, ni = 0; ((j < nX) && (ni < nnSize)); j++) {
			if (j != i) {
				double Lij = std::pow(X[i *mX +j], 2);
				if (Lij <= Li) {
					unsigned int ij = zi *nnSize +ni;
					P[ij] += std::exp(-Bi *Lij) /Ci;
					P[ij] += std::exp(-B[j *3 +0] *Lij) /B[j *3 +1];
					zP += P[ij];
					W[ij] = j;
					ni++;
				}
			}
		}
	}
}

// transform input euclidean-distances into probabilities
// FROM SPARSE-MATRIX DATA
void affMtx::mt_S2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W)
{
	for (unsigned int i = z_ini, zi = 0; i < z_end; i++, zi++) {
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[i *mX +v];
		double Bi = B[i *3 +0];
		double Ci = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int j = 0, ni = 0; ((j < nX) && (ni < nnSize)); j++) {
			if (j != i) {
				double Xj[mX];
				for(unsigned int v = 0; v < mX; v++) Xj[v] = X[j *mX +v];
				double Lij = spDist(mX, Xi, Xj);
				if (Lij <= Li) {
					unsigned int ij = zi *nnSize +ni;
					P[ij] += std::exp(-Bi *Lij) /Ci;
					P[ij] += std::exp(-B[j *3 +0] *Lij) /B[j *3 +1];
					zP += P[ij];
					W[ij] = j;
					ni++;
				}
			}
		}
	}
}
