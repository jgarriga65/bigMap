/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
//#include <algorithm>

#include "sqdist.h"
#include "affmtx.h"
#include "qtree.h"
#include "tsne.h"

// using namespace arma;
using namespace Rcpp;
using namespace std;

// DBL_EPSILON = 1.0e-9 (DBL_MIN = 1.0e-37)

// SHOOT OFF CONDITION CONTROL
// . minL1
// static const double minAtrForce = 1.0e-6 /(1.0 +2.0e-12);
// .minL2 (this one makes sense but does not work either)
static const double minAtrForce = std::sqrt(DBL_EPSILON) /(1.0 +2.0 *DBL_EPSILON);
// .minL3
// static const double minAtrForce = 2.0 *std::sqrt(DBL_EPSILON) /(1.0 +2.0 *DBL_EPSILON);

// thread affinity matrix
// [[Rcpp::export]]
double thread_affMtx(unsigned int z_ini, unsigned int z_end, SEXP sexpX, bool isDistance, bool isSparse, SEXP sexpB, SEXP sexpP, SEXP sexpW)
{
	// thread-size
	unsigned int z = z_end -z_ini;
	// recasting
	index_type offset;
	// P matrix
	BigMatrix *bigP = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpP));
	offset = bigP ->nrow() *bigP ->col_offset();
	double* thread_P = reinterpret_cast <double*> (bigP ->matrix()) +offset;
	unsigned int nnSize = bigP ->nrow();
	// W matrix
	BigMatrix *bigW = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpW));
	offset = bigW ->nrow() *bigW ->col_offset();
	unsigned int* thread_W = reinterpret_cast <unsigned int*> (bigW ->matrix()) +offset;
	// thread indexes
	int* zIdx = (int*) malloc(1 *sizeof(int));
	//
	affMtx* affmtx = new affMtx(sexpX, sexpB, zIdx, z, nnSize);
	if (isDistance)
		affmtx ->bh_D2P(z_ini, z_end, thread_P, thread_W);
	else if (isSparse)
		affmtx ->bh_S2P(z_ini, z_end, thread_P, thread_W);
	else
		affmtx ->bh_X2P(z_ini, z_end, thread_P, thread_W);
	// P normalization factor
	double zP = affmtx ->zP;
	// free memory
	delete affmtx;
	free(zIdx); zIdx = NULL;
	bigP = NULL; thread_P = NULL;
	bigW = NULL; thread_W = NULL;
	//
	return zP;
}

// compute thread repF and Q normalization factor
// [[Rcpp::export]]
double thread_repF(unsigned int z_ini, unsigned int z_end, SEXP sexpY, double theta, SEXP sexpR)
{
	// thread-size
	unsigned int z = z_end -z_ini;
	// current embedding
	Rcpp::XPtr<BigMatrix> bigY(sexpY);
	MatrixAccessor<double> mtxY(*bigY);
	unsigned int nX = bigY ->nrow();
	unsigned int mY = bigY ->ncol();
	// make local copy
	double* Y = (double*) malloc(nX *mY *sizeof(double));
	for (unsigned int i = 0, k = 0; i < nX; i++) {
		for(unsigned int d = 0; d < mY; d++, k++) Y[k] = mtxY[d][i];
	}
	// repulsive forces
	BigMatrix *bigR = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpR));
	index_type offset = bigR ->nrow() *bigR ->col_offset();
	double* repF = reinterpret_cast <double*> (bigR ->matrix()) +offset;
	for (unsigned int k = 0; k < z *mY; k++) repF[k] = .0;
	//
	double zQ = .0;
	if (theta == .0) {
		for (unsigned int i = z_ini, zi = 0; i < z_end; i++, zi++) {
			for (unsigned int j = 0; j < nX; j++) {
				double L[mY];
				double Lij = 1.0;
				for(unsigned int d = 0, ki = i *mY, kj = j *mY; d < mY; d++, ki++, kj++) {
					L[d] = Y[ki] -Y[kj];
					Lij += L[d] *L[d];
				}
				double Qij = 1.0 /Lij;
				for(unsigned int d = 0; d < mY; d++) repF[zi *mY +d] += Qij *Qij *L[d];
				zQ += Qij;
			}
		}
	} else {
		Quadtree* qtree = new Quadtree(Y, nX, mY, theta);
		for (unsigned int zi = z_ini, i = 0; zi < z_end; zi++, i++) {
			qtree ->repF(Y, &repF[i *mY], &zQ, zi);
		}
		delete(qtree);
	}
	// free memory
	free(Y); Y = NULL;
	return zQ;
}

// mt-SNE
// [[Rcpp::export]]
Rcpp::List thread_mIter(unsigned int z_ini, unsigned int z_end, SEXP sexpP, SEXP sexpW, SEXP sexpY, double sumP, double sumQ, SEXP sexpR, SEXP sexpU, SEXP sexpG, double eta, double alpha, double gain)
{
	// thread-size
	unsigned int z = z_end -z_ini;
	// recasting
	index_type offset;
	// P matrix
	BigMatrix *bigP = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpP));
	offset = bigP ->nrow() *bigP ->col_offset();
	double* P = reinterpret_cast <double*> (bigP ->matrix()) +offset;
	unsigned int nnSize = bigP ->nrow();
	// W matrix
	BigMatrix *bigW = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpW));
	offset = bigW ->nrow() *bigW ->col_offset();
	unsigned int* W = reinterpret_cast <unsigned int*> (bigW ->matrix()) +offset;
	// current embedding
	Rcpp::XPtr<BigMatrix> bigY(sexpY);
	MatrixAccessor<double> mtxY(*bigY);
	unsigned int nX = bigY ->nrow();
	unsigned int mY = bigY ->ncol();
	// make local copy
	double* Y = (double*) malloc(nX *mY *sizeof(double));
	for (unsigned int i = 0, k = 0; i < nX; i++) {
		for(unsigned int d = 0; d < mY; d++, k++) {
			Y[k] = mtxY[d][i];
		}
	}
	// repulsive forces
	BigMatrix *bigR = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpR));
	offset = bigR ->nrow() *bigR ->col_offset();
	double* repF = reinterpret_cast <double*> (bigR ->matrix()) +offset;
	// point-wise updates
	BigMatrix *bigU = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpU));
	offset = bigU ->nrow() *bigU ->col_offset();
	double* U = reinterpret_cast <double*> (bigU ->matrix()) +offset;
	// point-wise gains
	BigMatrix *bigG = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpG));
	offset = bigG ->nrow() *bigG ->col_offset();
	double* G = reinterpret_cast <double*> (bigG ->matrix()) +offset;
	// updated mapping positions
	Rcpp::NumericVector new_Y (z *mY);
	for (unsigned int k = 0, ij = z_ini *mY; k < z *mY; k++, ij++) new_Y[k] = Y[ij];
	// +++ TSNE
	double zCost = .0;
	for (unsigned int i = z_ini, zi = 0; i < z_end; i++, zi++) {
		double atrF[mY];
		for (unsigned int d = 0; d < mY; d++) atrF[d] = .0;
		// attractive force
		for (unsigned int ni = 0; ni < nnSize; ni++) {
			unsigned int j = W[zi *nnSize +ni];
			double L[mY];
			double Lij = 1.0;
			for(unsigned int d = 0, ki = i *mY, kj = j *mY; d < mY; d++, ki++, kj++) {
				L[d] = Y[ki] -Y[kj];
				Lij += L[d] *L[d];
			}
			unsigned int ij = zi *nnSize +ni;
			for(unsigned int d = 0; d < mY; d++) atrF[d] += P[ij] *L[d] /Lij;
			zCost -= P[ij] *std::log(Lij);
		}
		// update position
		for (unsigned int d = 0, k = zi *mY; d < mY; d++, k++) {
			// Att! atrF is double atrF[mY] !!!
			// if (std::abs(atrF[d] /sumP) > minAtrForce) {
				double dY = atrF[d] /sumP -repF[k] /sumQ;
				if (gain > .0) {
					if (signbit(dY) != signbit(U[k])) G[k] += (gain *.1);
					else G[k] *= (1.0 -gain /10.0);
				}
				U[k] = alpha *U[k] -eta *G[k] *dY;
				new_Y[k] += U[k];
			// }
		}
	}
	// free memory
	free(Y); Y = NULL;
	P = NULL;
	W = NULL;
	U = NULL;
	G = NULL;
	repF = NULL;
	//
	return Rcpp::List::create(Named("newY") = new_Y, Named("zCost") = zCost);
}
