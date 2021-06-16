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

#include "affmtx.h"
#include "tsne.h"

using namespace Rcpp;

// +++ SCKT -------------------------------------------------------------------

// Exact/Approx. SOCKET implementation of t-SNE
// [[Rcpp::export]]
double sckt_zTSNE(unsigned int thread_rank, unsigned int threads, unsigned int layers, SEXP sexpX, SEXP sexpB, SEXP sexpY, SEXP sexpI, arma::Col<double> &eRange, int iters, double nnSize, double lRate, double theta, double alpha, bool isDistance, bool isSparse, double latEx)
{
	// current mapping positions
	Rcpp::XPtr<BigMatrix> bigY(sexpY);
	MatrixAccessor<double> Y(*bigY);
	// sampled row indexes
	Rcpp::XPtr<BigMatrix> bigI(sexpI);
	MatrixAccessor<int> I(*bigI);
	unsigned int nX = bigI -> nrow();

	// get chunk breaks
	std::vector<int> breaks (threads +1, 0);
	for (unsigned int b = 0; b < threads; b++) breaks[b] = (int) b *(nX +1.0) /threads;
	breaks[threads] = nX;
	// get thread-size
	unsigned int thread_size = 0;
	for (unsigned int l = 0; l < layers; l++) {
		unsigned int a = (thread_rank + l) % threads;
		thread_size += (breaks[a+1] - breaks[a]);
	}

	// get thread-indexes & thread-layers
	int* zIdx = (int*) malloc(thread_size *sizeof(int));
	int* zLay = (int*) malloc(thread_size *sizeof(int));
	for (unsigned int l = 0, k = 0; l < layers; l++) {
		unsigned int a = (thread_rank + l) % threads;
		for (unsigned int i = breaks[a]; i < breaks[a+1]; i++, k++) {
			zIdx[k] = I[0][i];
			zLay[k] = l;
		}
	}

	// . cost function value
	double thread_Cost = 1.0;

	// . input affinities (P distribution)
	double* thread_P = (double*) calloc(thread_size *(thread_size -1) /2, sizeof(double));
	// . indexes of data-point pairs with significant atractive forces
	int* thread_W = (int*) calloc(thread_size *(thread_size -1) /2, sizeof(int));
	// . exaggeration factor
	double* thread_E = (double*) malloc(thread_size *sizeof(double));

	// +++ Compute affinity matrix
	affMtx* affmtx = new affMtx(sexpX, sexpB, zIdx, thread_size, nnSize, latEx);
	affmtx ->exgg(thread_E);
	if (isDistance)
		affmtx ->D2P(thread_P, thread_W);
	else if (isSparse)
		affmtx ->S2P(thread_P, thread_W);
	else
		affmtx ->X2P(thread_P, thread_W);

	// . starting mapping positions
	double* thread_Y = (double*) malloc(thread_size *2 *sizeof(double));
	for (unsigned int i = 0, k = 0; i < thread_size; i++) {
		for (unsigned int d = 0; d < 2; d++, k++){
			thread_Y[k] = Y[zLay[i] *2 +d][zIdx[i]];
		}
	}

	// +++ TSNE instance
	TSNE* tsne = new TSNE(thread_size, affmtx ->w, 2, eRange.begin(), iters, lRate, theta, alpha, affmtx ->zP, thread_E);
	// +++ run t-SNE
	tsne ->run2D(thread_P, thread_W, thread_Y);
	// +++ compute cost
	thread_Cost = tsne ->getCost(thread_P, thread_Y);

	// update mapping positions
	for (unsigned int l = 0, k = 0; l < layers; l++) {
		for (unsigned int a = (thread_rank + l) % threads, b = 0; b < (breaks[a+1] -breaks[a]); b++) {
			unsigned int j = I[0][breaks[a] +b];
			for (unsigned int d = 0; d < 2; d++, k++) Y[l *2 +d][j] = thread_Y[k];
		}
	}

	// free memory
	delete tsne;
	free(thread_Y); thread_Y = NULL;
	//
	delete affmtx;
	free(thread_E); thread_E = NULL;
	free(thread_W); thread_W = NULL;
	free(thread_P); thread_P = NULL;
	free(zIdx); zIdx = NULL;
	free(zLay); zLay = NULL;
	//
	return thread_Cost;
}

// +++ MPI --------------------------------------------------------------------

// get chunks by thread: indexes and mapping-positions
// [[Rcpp::export]]
void zChnks(Rcpp::List& Z_list, const arma::Mat<double>& Y, const arma::Col<int>& I, const Rcpp::List& brks_list)
{
	for (unsigned int z = 0; z < brks_list.length(); z++)
	{
		arma::Mat<int> brks = brks_list[z];
		arma::Mat<double> zChnk = Z_list[z];
		for (unsigned int l = 0, k = 0; l < brks.n_rows; l++) {
			for (unsigned int j1 = l *2, j2 = l *2 +1, i = brks(l, 0); i < brks(l, 1); i++, k++) {
				zChnk(k, 0) = I[i];
				zChnk(k, 1) = Y(I[i], j1);
				zChnk(k, 2) = Y(I[i], j2);
			}
		}
		Z_list[z] = zChnk;
	}
}

// restructure global mapping
// [[Rcpp::export]]
void updateY(arma::Mat<double>& Y, const arma::Col<int>& I, const Rcpp::List& zMap_list, const Rcpp::List& brks_list)
{
	for (unsigned int z = 0; z < zMap_list.length(); z++) {
		arma::Mat<int> thrd_brks = brks_list[z];
		arma::Mat<double> zY = zMap_list[z];
		for (unsigned int l = 0, k = 0; l < thrd_brks.n_rows; l++){
			for (unsigned int j1 = l *2, j2 = l *2 +1, i = thrd_brks(l, 0); i < thrd_brks(l, 1); i++, k++) {
				Y(I[i], j1) = zY(k, 0);
				Y(I[i], j2) = zY(k, 1);
			}
		}
	}
}

// Exact/Approx. MPI implementation of t-SNE
// [[Rcpp::export]]
double mpi_zTSNE(SEXP sexpX, SEXP sexpB, arma::Mat<double> &Y, arma::Col<int> indexes, arma::Col<double> &eRange, int iters, double nnSize, double lRate, double theta, double alpha, bool isDistance, bool isSparse, double latEx)
{
	// unsigned int N = bmL->nrow();
	unsigned int thread_size = Y.n_rows;
	// int* zIdx = reinterpret_cast <int*> (I.begin());
	int* zIdx = indexes.begin();
	// . cost function value
	double thread_Cost = 1.0;

	// . input affinities (P distribution)
	double* thread_P = (double*) calloc(thread_size * (thread_size -1) /2, sizeof(double));
	// . indexes of data-point pairs with significant atractive forces
	int* thread_W = (int*) calloc(thread_size * (thread_size -1) /2, sizeof(int));
	// . exaggeration factor
	double* thread_E = (double*) malloc(thread_size *sizeof(double));

	// +++ Compute affinity matrix
	affMtx* affmtx = new affMtx(sexpX, sexpB, zIdx, thread_size, nnSize, latEx);
	affmtx ->exgg(thread_E);
	if (isDistance)
		affmtx ->D2P(thread_P, thread_W);
	else if (isSparse)
		affmtx ->S2P(thread_P, thread_W);
	else
		affmtx ->X2P(thread_P, thread_W);

	// . starting mapping positions
	double* thread_Y = (double*) malloc(thread_size *2 *sizeof(double));
	for (unsigned int i = 0, k = 0; i < thread_size; i++) {
		for (unsigned int d = 0; d < 2; d++, k++) thread_Y[k] = Y(i, d);
	}

	// +++ TSNE instance
	TSNE* tsne = new TSNE(thread_size, affmtx ->w, 2, eRange.begin(), iters, lRate, theta, alpha, affmtx ->zP, thread_E);
	// +++ run t-SNE
	tsne ->run2D(thread_P, thread_W, thread_Y);
	// +++ compute cost
	thread_Cost = tsne ->getCost(thread_P, thread_Y);

	// update mapping positions
	for (unsigned int i = 0, k = 0; i < thread_size; i++) {
		for (unsigned int d = 0; d < 2; d++, k++) Y(i, d) = thread_Y[k];
	}

	// free memory
	delete tsne;
	free(thread_Y); thread_Y = NULL;
	// deallocate memory
	delete affmtx;
	free(thread_E); thread_E = NULL;
	free(thread_W); thread_W = NULL;
	free(thread_P); thread_P = NULL;
	zIdx = NULL;
	//
	return thread_Cost;
}
