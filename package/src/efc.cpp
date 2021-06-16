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
	int* thread_W = reinterpret_cast <int*> (bigW ->matrix()) +offset;
	// thread indexes
	int* zIdx = (int*) malloc(1 *sizeof(int));
	//
	affMtx* affmtx = new affMtx(sexpX, sexpB, zIdx, z, nnSize, 1.0);
	if (isDistance)
		affmtx ->efc_D2P(z_ini, z_end, thread_P, thread_W);
	else if (isSparse)
		affmtx ->efc_S2P(z_ini, z_end, thread_P, thread_W);
	else
		affmtx ->efc_X2P(z_ini, z_end, thread_P, thread_W);
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

// compute thread Q normalization factor
// [[Rcpp::export]]
double thread_qNorm(unsigned int z_ini, unsigned int z_end, SEXP sexpY, double theta)
{
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
	//
	double zQ = .0;
	if (theta == .0) {
		for (unsigned int zi = z_ini; zi < z_end; zi++) {
			for (unsigned int zj = 0; zj < nX; zj++) {
				double Lij = 1.0;
				for (unsigned int d = 0; d < mY; d++) Lij += std::pow(Y[zi *mY +d] -Y[zj *mY +d], 2.0);
				zQ += 1.0 /Lij;
			}
		}
	} else {
		double repF[mY];
		Quadtree* qtree = new Quadtree(Y, nX, mY);
		for (unsigned int zi = z_ini; zi < z_end; zi++) qtree ->repF(Y, repF, theta, &zQ, zi);
		delete(qtree);
	}
	// free memory
	free(Y); Y = NULL;
	return zQ;
}

// embedding final compression
// [[Rcpp::export]]
Rcpp::NumericVector thread_itrRun(unsigned int z_ini, unsigned int z_end, SEXP sexpP, SEXP sexpW, SEXP sexpY, double zP, double zQ, double lRate, double theta)
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
	int* W = reinterpret_cast <int*> (bigW ->matrix()) +offset;
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
	// updated mapping positions
	Rcpp::NumericVector new_Y (z *mY);
	for (unsigned int k = 0, ij = z_ini *mY; k < z *mY; k++, ij++) new_Y[k] = Y[ij];
	// learning-rate
	double eta = lRate *4.0;
	// +++ TSNE
	if (theta == .0) {
		for (unsigned int zi = z_ini, i = 0; zi < z_end; zi++, i++){
			double atrF[mY];
			double repF[mY];
			for (unsigned int d = 0; d < mY; d++) {
				atrF[d] = .0;
				repF[d] = .0;
			}
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
			}
			// repulsive force
			for (unsigned int zj = 0; zj < nX; zj++) {
				double L[mY];
				double Lij = 1.0;
				for(unsigned int d = 0; d < mY; d++) {
					L[d] = Y[zi *mY +d] -Y[zj *mY +d];
					Lij += L[d] *L[d];
				}
				double Qij = 1.0 /Lij;
				for(unsigned int d = 0; d < mY; d++) repF[d] += Qij *Qij *L[d];
			}
			// update position
			for (unsigned int d = 0; d < mY; d++){
				new_Y[i *mY +d] -= eta *(atrF[d] /zP -repF[d] /zQ);
			}
		}
	} else {
		double repF[mY];
		double zQi = .0;
		Quadtree* qtree = new Quadtree(Y, nX, mY);
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
				for(unsigned int d = 0; d < mY; d++) {
					if (P[ij] > 0) {
						atrF[d] += P[ij] *Qij *L[d];
					}
				}
			}
			// repulsive force
			qtree->repF(Y, repF, theta, &zQi, zi);
			// update position
			// new_Y[i *mY +0] = atrF[0] /zP;
			// new_Y[i *mY +1] = repF[0] /zQ;
			for (unsigned int d = 0; d < mY; d++){
				new_Y[i *mY +d] -= eta *(atrF[d] /zP -repF[d] /zQ);
			}
		}
		delete qtree;
	}
	// free memory
	free(Y); Y = NULL;
	bigP = NULL; P = NULL;
	bigW = NULL; W = NULL;
	//
	return new_Y;
}
