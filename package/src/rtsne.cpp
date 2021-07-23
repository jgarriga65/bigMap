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
		affmtx ->efr_D2P(z_ini, z_end, thread_P, thread_W);
	else if (isSparse)
		affmtx ->efr_S2P(z_ini, z_end, thread_P, thread_W);
	else
		affmtx ->efr_X2P(z_ini, z_end, thread_P, thread_W);
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
	double minL = DBL_EPSILON /4.0;
	//
	double zQ = .0;
	if (theta == .0) {
		for (unsigned int zi = z_ini, i = 0; zi < z_end; zi++, i++) {
			for (unsigned int zj = 0; zj < nX; zj++) {
				double L[mY];
				double Lij = 1.0;
				for(unsigned int d = 0; d < mY; d++) {
					L[d] = Y[zi *mY +d] -Y[zj *mY +d];
					if (std::abs(L[d]) < minL) {
						L[d] = (Y[zi *mY +d] >Y[zj *mY +d]) ? minL : -minL;
					}
					Lij += L[d] *L[d];
				}
				double Qij = 1.0 /Lij;
				for(unsigned int d = 0; d < mY; d++) repF[i *d] += Qij *Qij *L[d];
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

// embedding final compression
// [[Rcpp::export]]
Rcpp::List thread_iter(unsigned int z_ini, unsigned int z_end, SEXP sexpP, SEXP sexpW, SEXP sexpY, double sumP, double sumQ, SEXP sexpR, SEXP sexpU, SEXP sexpG, double lRate, double alpha)
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
		for(unsigned int d = 0; d < mY; d++, k++) Y[k] = mtxY[d][i];
	}
	// point-wise updates
	BigMatrix *bigU = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpU));
	offset = bigU ->nrow() *bigU ->col_offset();
	double* U = reinterpret_cast <double*> (bigU ->matrix()) +offset;
	// point-wise gains
	BigMatrix *bigG = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpG));
	offset = bigG ->nrow() *bigG ->col_offset();
	double* G = reinterpret_cast <double*> (bigG ->matrix()) +offset;
	// repulsive forces
	BigMatrix *bigR = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpR));
	offset = bigR ->nrow() *bigR ->col_offset();
	double* repF = reinterpret_cast <double*> (bigR ->matrix()) +offset;
	// updated mapping positions
	Rcpp::NumericVector new_Y (z *mY);
	for (unsigned int k = 0, ij = z_ini *mY; k < z *mY; k++, ij++) new_Y[k] = Y[ij];
	// learning-rate
	double eta = lRate *4.0;
	double zCost = .0;
	// +++ TSNE
	for (unsigned int zi = z_ini, i = 0; zi < z_end; zi++, i++){
		double atrF[mY];
		for (unsigned int d = 0; d < mY; d++) atrF[d] = .0;
		// attractive force
		for (unsigned int ni = 0; ni < nnSize; ni++) {
			unsigned int zj = W[i *nnSize +ni];
			double L[mY];
			double Lij = 1.0;
			for(unsigned int d = 0; d < mY; d++) {
				L[d] = Y[zi *mY +d] -Y[zj *mY +d];
				Lij += L[d] *L[d];
			}
			double Qij = 1.0 /Lij;
			unsigned int ij = i *nnSize +ni;
			for(unsigned int d = 0; d < mY; d++) atrF[d] += P[ij] *Qij *L[d];
			zCost -= P[ij] *std::log(Lij);
		}
		// update position
		// new_Y[i *mY + 0] = atrF[0] /sumP;
		// new_Y[i *mY + 1] = repF[i *mY +0] /sumQ;
		for (unsigned int d = 0; d < mY; d++){
			double dY = atrF[d] /sumP -repF[i *mY +d] /sumQ;
			if (signbit(dY) != signbit(U[i *mY +d])) {
				G[i *mY +d] += .2;
			}
			else {
				G[i *mY +d] *= .8;
				G[i *mY +d] = std::max(G[i *mY +d], .01);
			}
			U[i *mY +d] = alpha *U[i *mY + d] -eta *G[i *mY +d] *dY;
			new_Y[i *mY +d] += U[i *mY +d];
		}
	}
	// free memory
	free(Y); Y = NULL;
	P = NULL;
	W = NULL;
	U = NULL;
	repF = NULL;
	//
	return Rcpp::List::create(Named("newY") = new_Y, Named("zCost") = zCost);
}
