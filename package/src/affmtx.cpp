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
affMtx::affMtx(SEXP sexpX, SEXP sexpB, int* zIdx, unsigned int z, unsigned int nnSize, double latEx) : zIdx(zIdx), z(z), nnSize(nnSize), latEx(latEx)
{

	// Att!! X and B transposed
	index_type offset;

	BigMatrix *bigX = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpX));
	offset = bigX ->nrow() *bigX ->col_offset();
	this ->X = reinterpret_cast <double*> (bigX ->matrix()) +offset;
	this ->nX = bigX ->ncol();
	this ->mX = bigX ->nrow();
	this ->w = 0;

	BigMatrix *bigB = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpB));
	offset = bigB ->nrow() *bigB ->col_offset();
	this ->B = reinterpret_cast <double*> (bigB ->matrix()) +offset;

	// thread fraction of affinities
	this ->zP = .0;

	// clean
	bigX = NULL;
	bigB = NULL;
}

// local exaggeration factor
void affMtx::exgg(double* E)
{
	for (unsigned int i = 0; i < z; i++) E[i] = B[zIdx[i] *4 +3];
}

// transform input similarities into probabilities
// FROM INPUT-DATA
void affMtx::X2P(double* P, int* W)
{
	std::vector<int> nn(z, 0);
	for (unsigned int i = 0, ij = 0; i < z; i++) {
		unsigned int zi = zIdx[i];
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[zi *mX +v];
		double Bi = B[zi *4 +0];
		double Zi = B[zi *4 +1];
		double Li = B[zi *4 +2];
		for (unsigned int j = i +1; j < z; j++, ij++) {
			unsigned int zj = zIdx[j];
			double Lij = .0;
			for(unsigned int v = 0; v < mX; v++) Lij += pow(Xi[v] -X[zj *mX +v], 2);
			double Bj = B[zj *4 +0];
			double Zj = B[zj *4 +1];
			double Lj = B[zj *4 +2];
			if ((Lij <= Li) || (Lij <= Lj)) {
				if ((nn[i] < nnSize) || (nn[j] < nnSize)) {
					P[ij] += std::exp(-Bi *Lij) /Zi;
					P[ij] += std::exp(-Bj *Lij) /Zj;
					zP += P[ij];
					nn[i] ++;
					nn[j] ++;
					W[w] = i *z +j;
					w ++;
				}
			}
		}
	}
	zP *= 2.0 /latEx;
}

// transform input similarities into probabilities
// FROM FULL-DISTANCE-MATRIX
void affMtx::D2P(double* P, int* W)
{
	std::vector<int> nn(z, 0);
	for (unsigned int i = 0, ij = 0; i < z; i++) {
		unsigned int zi = zIdx[i];
		double Bi = B[zi *4 +0];
		double Zi = B[zi *4 +1];
		double Li = B[zi *4 +2];
		for (unsigned int j = i +1; j < z; j++, ij++) {
			unsigned int zj = zIdx[j];
			double Lij = std::pow(X[zi *mX +zj], 2);
			double Bj = B[zj *4 +0];
			double Zj = B[zj *4 +1];
			double Lj = B[zj *4 +2];
			if ((Lij <= Li) || (Lij <= Lj)) {
				if ((nn[i] < nnSize) || (nn[j] < nnSize)) {
					P[ij] += std::exp(-Bi *Lij) /Zi;
					P[ij] += std::exp(-Bj *Lij) /Zj;
					zP += P[ij];
					nn[i] ++;
					nn[j] ++;
					W[w] = i *z +j;
					w ++;
				}
			}
		}
	}
	zP *= 2.0 /latEx;
}

// transform input euclidean-distances into probabilities
// FROM SPARSE-MATRIX DATA
void affMtx::S2P(double* P, int* W)
{
	std::vector<int> nn(z, 0);
	for (unsigned int i = 0, ij = 0; i < z; i++) {
		unsigned int zi = zIdx[i];
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[zi *mX +v];
		double Bi = B[zi *4 +0];
		double Zi = B[zi *4 +1];
		double Li = B[zi *4 +2];
		for (unsigned int j = i +1; j < z; j++, ij++) {
			unsigned int zj = zIdx[j];
			double Xj[mX];
			for(unsigned int v = 0; v < mX; v++) Xj[v] = X[zj *mX +v];
			double Lij = spDist(mX, Xi, Xj);
			double Bj = B[zj *4 +0];
			double Zj = B[zj *4 +1];
			double Lj = B[zj *4 +2];
			if ((Lij <= Li) || (Lij <= Lj)) {
				if ((nn[i] < nnSize) || (nn[j] < nnSize)) {
					P[ij] += std::exp(-Bi *Lij) /Zi;
					P[ij] += std::exp(-Bj *Lij) /Zj;
					zP += P[ij];
					nn[i] ++;
					nn[j] ++;
					W[w] = i *z +j;
					w ++;
				}
			}
		}
	}
	zP *= 2.0 /latEx;
}

// +++++++++++++++++++++++++++++++++++ embedding final compression

// transform input similarities into probabilities
// FROM INPUT-DATA
void affMtx::efr_X2P(unsigned int z_ini, unsigned int z_end, double* P, int* W)
{
	for (unsigned int zi = z_ini, i = 0; zi < z_end; zi++, i++) {
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[zi *mX +v];
		double Bi = B[zi *4 +0];
		double Zi = B[zi *4 +1];
		double Li = B[zi *4 +2];
		for (unsigned int zj = 0, ni = 0; ((zj < nX) && (ni < nnSize)); zj++) {
			if (zj != zi) {
				double Lij = .0;
				for(unsigned int v = 0; v < mX; v++) Lij += pow(Xi[v] -X[zj *mX +v], 2);
				if (Lij <= Li) {
					double Bj = B[zj *4 +0];
					double Zj = B[zj *4 +1];
					unsigned int ij = i *nnSize +ni;
					P[ij] += std::exp(-Bi *Lij) /Zi;
					P[ij] += std::exp(-Bj *Lij) /Zj;
					zP += P[ij];
					W[ij] = zj;
					ni ++;
				}
			}
		}
	}
	//zP *= 2.0;	// this is not triangular !!!
}


// transform input similarities into probabilities
// FROM FULL-DISTANCE-MATRIX
void affMtx::efr_D2P(unsigned int z_ini, unsigned int z_end, double* P, int* W)
{
	for (unsigned int zi = z_ini, i = 0; zi < z_end; zi++, i++) {
		double Bi = B[zi *4 +0];
		double Zi = B[zi *4 +1];
		double Li = B[zi *4 +2];
		for (unsigned int zj = 0, ni = 0; ((zj < nX) && (ni < nnSize)); zj++) {
			if (zj != zi) {
				double Lij = std::pow(X[zi *mX +zj], 2);
				if (Lij <= Li) {
					double Bj = B[zj *4 +0];
					double Zj = B[zj *4 +1];
					unsigned int ij = i *nnSize +ni;
					P[ij] += std::exp(-Bi *Lij) /Zi;
					P[ij] += std::exp(-Bj *Lij) /Zj;
					zP += P[ij];
					W[ij] = zj;
					ni ++;
				}
			}
		}
	}
}

// transform input euclidean-distances into probabilities
// FROM SPARSE-MATRIX DATA
void affMtx::efr_S2P(unsigned int z_ini, unsigned int z_end, double* P, int* W)
{
	for (unsigned int zi = z_ini, i = 0; zi < z_end; zi++, i++) {
		double Xi[mX];
		for(unsigned int v = 0; v < mX; v++) Xi[v] = X[zi *mX +v];
		double Bi = B[zi *4 +0];
		double Zi = B[zi *4 +1];
		double Li = B[zi *4 +2];
		for (unsigned int zj = 0, ni = 0; ((zj < nX) && (ni < nnSize)); zj++) {
			double Xj[mX];
			for(unsigned int v = 0; v < mX; v++) Xj[v] = X[zj *mX +v];
			double Lij = spDist(mX, Xi, Xj);
			if (Lij <= Li) {
				double Bj = B[zj *4 +0];
				double Zj = B[zj *4 +1];
				unsigned int ij = i *nnSize +ni;
				P[ij] += std::exp(-Bi *Lij) /Zi;
				P[ij] += std::exp(-Bj *Lij) /Zj;
				zP += P[ij];
				W[ij] = zj;
				ni ++;
			}
		}
	}
}
