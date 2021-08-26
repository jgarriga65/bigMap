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

	// thread fraction of affinities
	this ->zP = .0;

	// clean
	bigX = NULL;
	bigB = NULL;
}

// transform input similarities into probabilities
// FROM INPUT-DATA
void affMtx::X2P(double* P, unsigned int* W)
{
	std::vector<unsigned int> n(z, 0);
	for (unsigned int zi = 0; zi < z; zi++) {
		unsigned int i = zIdx[zi];
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[i *mX +v];
		double Bi = B[i *3 +0];
		double Zi = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int zj = zi +1; zj < z; zj++) {
			unsigned int j = zIdx[zj];
			double Lij = .0;
			for(unsigned int v = 0; v < mX; v++) Lij += pow(Xi[v] -X[j *mX +v], 2);
			double Bj = B[j *3 +0];
			double Zj = B[j *3 +1];
			double Lj = B[j *3 +2];
			if ((Lij <= Li) && (n[zi] < nnSize)) {
				unsigned int ij = zi *nnSize +n[zi];
				P[ij]  = std::exp(-Bi *Lij) /Zi;
				P[ij] += std::exp(-Bj *Lij) /Zj;
				zP += P[ij];
				W[ij] = zj;
				n[zi]++;
			}
			if ((Lij <= Lj) && (n[zj] < nnSize)) {
				unsigned int ji = zj *nnSize +n[zj];
				P[ji]  = std::exp(-Bi *Lij) /Zi;
				P[ji] += std::exp(-Bj *Lij) /Zj;
				zP += P[ji];
				W[ji] = zi;
				n[zj]++;
			}
		}
	}
}

// transform input similarities into probabilities
// FROM FULL-DISTANCE-MATRIX
void affMtx::D2P(double* P, unsigned int* W)
{
	std::vector<unsigned int> n(z, 0);
	for (unsigned int zi = 0; zi < z; zi++) {
		unsigned int i = zIdx[zi];
		double Bi = B[i *3 +0];
		double Zi = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int zj = zi +1; zj < z; zj++) {
			unsigned int j = zIdx[zj];
			double Lij = std::pow(X[i *mX +j], 2);
			double Bj = B[j *3 +0];
			double Zj = B[j *3 +1];
			double Lj = B[j *3 +2];
			if ((Lij <= Li) && (n[zi] < nnSize)) {
				unsigned int ij = zi *nnSize +n[zi];
				P[ij]  = std::exp(-Bi *Lij) /Zi;
				P[ij] += std::exp(-Bj *Lij) /Zj;
				zP += P[ij];
				W[ij] = zj;
				n[zi]++;
			}
			if ((Lij <= Lj) && (n[zj] < nnSize)) {
				unsigned int ji = zj *nnSize +n[zj];
				P[ji]  = std::exp(-Bi *Lij) /Zi;
				P[ji] += std::exp(-Bj *Lij) /Zj;
				zP += P[ji];
				W[ji] = zi;
				n[zj]++;
			}
		}
	}
}

// transform input euclidean-distances into probabilities
// FROM SPARSE-MATRIX DATA
void affMtx::S2P(double* P, unsigned int* W)
{
	std::vector<unsigned int> n(z, 0);
	for (unsigned int zi = 0; zi < z; zi++) {
		unsigned int i = zIdx[zi];
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[i *mX +v];
		double Bi = B[i *3 +0];
		double Zi = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int zj = zi +1; zj < z; zj++) {
			unsigned int j = zIdx[zj];
			double Xj[mX];
			for(unsigned int v = 0; v < mX; v++) Xj[v] = X[j *mX +v];
			double Lij = spDist(mX, Xi, Xj);
			double Bj = B[j *3 +0];
			double Zj = B[j *3 +1];
			double Lj = B[j *3 +2];
			if ((Lij <= Li) && (n[zi] < nnSize)) {
				unsigned int ij = zi *nnSize +n[zi];
				P[ij]  = std::exp(-Bi *Lij) /Zi;
				P[ij] += std::exp(-Bj *Lij) /Zj;
				zP += P[ij];
				W[ij] = zj;
				n[zi]++;
			}
			if ((Lij <= Lj) && (n[zj] < nnSize)) {
				unsigned int ji = zj *nnSize +n[zj];
				P[ji]  = std::exp(-Bi *Lij) /Zi;
				P[ji] += std::exp(-Bj *Lij) /Zj;
				zP += P[ji];
				W[ji] = zi;
				n[zj]++;
			}
		}
	}
}

// void affMtx::S2P_ordered(double* P, unsigned int* W)
// {
// 	// +++ possible take-home message
// 	// . when we have many neighbors at the same distance it can be beneficial
// 	// to take whatever information that might help to discern the most
// 	// appropiate set of neighbors like e.g. the ordinal sequence for PRIMES
//  	unsigned int r = 0;
// 	std::vector<int> rank(z);
// 	std::iota(rank.begin(), rank.end(), r++);
// 	std::sort(rank.begin(), rank.end(),
// 		[&](int a, int b) {return zIdx[a] < zIdx[b];}
// 	);
// 	//
// 	std::vector<unsigned int> n(z, 0);
// 	for (unsigned int ri = 0; ri < z; ri++) {
// 		unsigned int zi = rank[ri];
// 		unsigned int i = zIdx[zi];
// 		double Xi[mX];
// 		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[i *mX +v];
// 		double Bi = B[i *3 +0];
// 		double Zi = B[i *3 +1];
// 		double Li = B[i *3 +2];
// 		for (unsigned int rj = ri +1; rj < z; rj++) {
// 			unsigned int zj = rank[rj];
// 			unsigned int j = zIdx[zj];
// 			double Xj[mX];
// 			for(unsigned int v = 0; v < mX; v++) Xj[v] = X[j *mX +v];
// 			double Lij = spDist(mX, Xi, Xj);
// 			double Bj = B[j *3 +0];
// 			double Zj = B[j *3 +1];
// 			double Lj = B[j *3 +2];
// 			if ((Lij <= Li) && (n[zi] < nnSize)) {
// 				unsigned int ij = zi *nnSize +n[zi];
// 				P[ij]  = std::exp(-Bi *Lij) /Zi;
// 				P[ij] += std::exp(-Bj *Lij) /Zj;
// 				zP += P[ij];
// 				W[ij] = zj;
// 				n[zi]++;
// 			}
// 			if ((Lij <= Lj) && (n[zj] < nnSize)) {
// 				unsigned int ji = zj *nnSize +n[zj];
// 				P[ji]  = std::exp(-Bi *Lij) /Zi;
// 				P[ji] += std::exp(-Bj *Lij) /Zj;
// 				zP += P[ji];
// 				W[ji] = zi;
// 				n[zj]++;
// 			}
// 		}
// 	}
// }

// void affMtx::S2P_exact(double* P, unsigned int* W)
// {
// 	std::vector<double> Ki(z);
// 	for (unsigned int zi = 0; zi < z; zi++) {
// 		unsigned int i = zIdx[zi];
// 		double Xi[mX];
// 		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[i *mX +v];
// 		double Bi = B[i *3 +0];
// 		double Zi = B[i *3 +1];
// 		double Li = B[i *3 +2];
// 		for (unsigned int zj = 0; zj < z; zj++) {
// 			Ki[zj] = .0;
// 			unsigned int j = zIdx[zj];
// 			if (j != i) {
// 				double Xj[mX];
// 				for(unsigned int v = 0; v < mX; v++) Xj[v] = X[j *mX +v];
// 				double Lij = spDist(mX, Xi, Xj);
// 				if (Lij <= Li) Ki[zj] = std::exp(-Bi *Lij) /Zi +std::exp(-B[j *3] *Lij) /B[j *3 +1];
// 			}
// 		}
// 		row_affinity(z, nnSize, i, Ki.data(), P, W, zP);
// 	}
// }

// +++++++++++++++++++++++++++++++++++ embedding final compression

// transform input similarities into probabilities
// FROM INPUT-DATA
void affMtx::bh_X2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W)
{
	for (unsigned int i = z_ini, zi = 0; i < z_end; i++, zi++) {
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[i *mX +v];
		double Bi = B[i *3 +0];
		double Zi = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int j = 0, ni = 0; ((j < nX) && (ni < nnSize)); j++) {
			if (j != i) {
				double Lij = .0;
				for(unsigned int v = 0; v < mX; v++) Lij += std::pow(Xi[v] -X[j *mX +v], 2);
				if (Lij <= Li) {
					double Bj = B[j *3 +0];
					double Zj = B[j *3 +1];
					unsigned int ij = zi *nnSize +ni;
					P[ij] += std::exp(-Bi *Lij) /Zi;
					P[ij] += std::exp(-Bj *Lij) /Zj;
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
void affMtx::bh_D2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W)
{
	for (unsigned int i = z_ini, zi = 0; i < z_end; i++, zi++) {
		double Bi = B[i *3 +0];
		double Zi = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int j = 0, ni = 0; ((j < nX) && (ni < nnSize)); j++) {
			if (j != i) {
				double Lij = std::pow(X[i *mX +j], 2);
				if (Lij <= Li) {
					double Bj = B[j *3 +0];
					double Zj = B[j *3 +1];
					unsigned int ij = zi *nnSize +ni;
					P[ij] += std::exp(-Bi *Lij) /Zi;
					P[ij] += std::exp(-Bj *Lij) /Zj;
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
void affMtx::bh_S2P(unsigned int z_ini, unsigned int z_end, double* P, unsigned int* W)
{
	for (unsigned int i = z_ini, zi = 0; i < z_end; i++, zi++) {
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[i *mX +v];
		double Bi = B[i *3 +0];
		double Zi = B[i *3 +1];
		double Li = B[i *3 +2];
		for (unsigned int j = 0, ni = 0; ((j < nX) && (ni < nnSize)); j++) {
			if (j != i) {
				double Xj[mX];
				for(unsigned int v = 0; v < mX; v++) Xj[v] = X[j *mX +v];
				double Lij = spDist(mX, Xi, Xj);
				if (Lij <= Li) {
					double Bj = B[j *3 +0];
					double Zj = B[j *3 +1];
					unsigned int ij = zi *nnSize +ni;
					P[ij] += std::exp(-Bi *Lij) /Zi;
					P[ij] += std::exp(-Bj *Lij) /Zj;
					zP += P[ij];
					W[ij] = j;
					ni++;
				}
			}
		}
	}
}

// [[Rcpp::export]]
Rcpp::NumericVector sortZidx(unsigned int z, Rcpp::NumericVector zIdx)
{
 	unsigned int r = 0;
	std::vector<int> rank(z);
	std::iota(rank.begin(), rank.end(), r++);
	std::sort(rank.begin(), rank.end(),
		[&](int a, int b) {return zIdx[a] < zIdx[b];}
	);
	Rcpp::NumericVector zIdx_sorted(z);
	for (unsigned int r = 0; r < z; r++) zIdx_sorted[r] = zIdx[rank[r]];
	return zIdx_sorted;
}

/***R
check1 <- function(n, z)
{
	zIdx <- sample(seq(n))[1:z]
	zIdx_sorted <- sortZidx(z, zIdx)
	zIdx_ret <- cbind(zIdx, zIdx_sorted)
	return(zIdx_ret)
}

*/
