/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>

#include <numeric>      // std::iota
#include <algorithm>    // std::nth_element, std::sort, std::stable_sort

#include "sqdist.h"
#include "affmtx.h"

// using namespace arma;
using namespace Rcpp;
using namespace std;

// check neighbors' set size
// [[Rcpp::export]]
Rcpp::NumericVector nnSS_chk(SEXP sexpX, SEXP sexpB, arma::Col<int> indexes, bool isDistance, bool isSparse, unsigned int nnSize)
{
	// indexes
	int* zIdx = indexes.begin();
	unsigned int thread_size = indexes.size();
	//. thread affinity matrix
	unsigned int aff_size = thread_size *(thread_size -1) /2;
	double* thread_P = (double*) calloc(aff_size, sizeof(double));
	// . indexes of data-point pairs with significant atractive forces
	unsigned int* thread_W = (unsigned int*) calloc(aff_size, sizeof(unsigned int));
	// . affinity matrix
	affMtx* affmtx = new affMtx(sexpX, sexpB, zIdx, thread_size, nnSize);
	if (isDistance)
		affmtx ->D2P(thread_P, thread_W);
	else if (isSparse)
		affmtx ->S2P(thread_P, thread_W);
	else
		affmtx ->X2P(thread_P, thread_W);
	// . neighbors' set size
	Rcpp::NumericVector thread_nnSS(thread_size);
	thread_nnSS.fill(0);
	for (unsigned int i = 0, ij = 0; i < thread_size; i++) {
		for (unsigned int j = i +1; j < thread_size; j++, ij++) {
			if (thread_P[ij] > 0) {
				thread_nnSS[i] ++;
				thread_nnSS[j] ++;
			}
		}
	}
	// free memory
	delete affmtx;
	free(thread_W); thread_W = NULL;
	free(thread_P); thread_P = NULL;
	zIdx = NULL;
	return thread_nnSS;
}
