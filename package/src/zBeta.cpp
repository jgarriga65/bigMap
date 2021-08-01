/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

// #include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <algorithm>

#include "sqdist.h"
#include "lambertW.h"

// using namespace arma;
using namespace Rcpp;
using namespace std;

// Center/Scale input-data matrix (avoids unclear numerical problems)
// Expects input-data as a transposed big.matrix !!
// [[Rcpp::export]]
void centerScale(SEXP sexpX, bool is_distance, bool is_sparse)
{
	BigMatrix *bigX = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpX));
	index_type offset = bigX ->nrow() *bigX ->col_offset();
	double* X = reinterpret_cast <double*> (bigX ->matrix()) +offset;
	unsigned int nX = bigX ->ncol();
	unsigned int mX = bigX ->nrow();
	// +++ Center
	if (!is_distance)
	{
		std::vector<double> xMean(mX, 0);
		// Att.!!! THIS MIGHT BE NECESSARY OUT OF THE PRIMES DATA SET
		// if (is_sparse)
		// {
		// 	std::vector<double> xSize(mX, 0);
		// 	for (unsigned int i = 0, ij = 1; i < nX; i++) {
		// 		for (unsigned int j = 1; j < mX; j++, j++, ij++, ij++) {
		// 			xMean[j] += X[ij];
		// 			if (X[ij] > 0) xSize[j] ++;
		// 		}
		// 	}
		// 	// for (unsigned int j = 1; j < mX; j++, j++) xMean[j] /= (double) nX;
		// 	for (unsigned int j = 1; j < mX; j++, j++) xMean[j] /= xSize[j];
		// 	// Subtract data mean
		// 	for (unsigned int i = 0, ij = 1; i < nX; i++) {
		// 		for (unsigned int j = 1; j < mX; j++, j++, ij++, ij++) X[ij] -= xMean[j];
		// 	}
		// }
		// else
		if (!is_sparse)
		{
			for (unsigned int i = 0, ij = 0; i < nX; i++) {
				for (unsigned int j = 0; j < mX; j++, ij++) xMean[j] += X[ij];
			}
			for (unsigned int j = 0; j < mX; j++) xMean[j] /= (double) nX;
			// Subtract data mean
			for (unsigned int i = 0, ij = 0; i < nX; i++) {
				for (unsigned int j = 0; j < mX; j++, ij++) X[ij] -= xMean[j];
			}
		}
	}
	// +++ Scale
	double Xmax = .0;
	if (is_sparse)
	{
		for (unsigned int ij = 1; ij < nX *mX; ij++, ij++) if (std::abs(X[ij]) > Xmax) Xmax = std::abs(X[ij]);
		for (unsigned int ij = 1; ij < nX *mX; ij++, ij++) X[ij] /= Xmax;
	}
	else
	{
		for (unsigned int ij = 0; ij < nX *mX; ij++) if (std::abs(X[ij]) > Xmax) Xmax = std::abs(X[ij]);
		for (unsigned int ij = 0; ij < nX *mX; ij++) X[ij] /= Xmax;
	}
	// +++ free memory
	bigX = NULL; X = NULL;
}

// entropy (up to qNN) (~ FIt-SNE version)
// Zi = 0 is not more a problem and Bi are almost identical to Fit-SNE (!!)
double entropy_1(std::vector<double> Li, double Bi, double &Zi, int qNN)
{
	// Att!! Hi = \sum_j Pij/Pi log(Pij/Pi); (no negative sign !!)
	// because:
	// Pij = exp(-Bi *Lij)
	// log(Pij) = -Bi *Lij (note the sign !!!)
	Zi = DBL_MIN;	// this setting along with (*) avoids Zi = 0;
	double Hi = .0;
	// +++ entropy (up to qNN)
	for (unsigned int j = 0; j < qNN; j++) {
		double Kij = std::exp(-Bi *Li[j]);
		if (Li[j] == 0) Kij = DBL_MIN;	// (*)
		Zi += Kij;
		Hi += Li[j] *Kij;
	}
	Hi *= Bi /Zi;
	Hi += std::log(Zi);
	return Hi;
}

// entropy (up to qNN) (my version)
double entropy(std::vector<double> Li, double Bi, double &Zi, int qNN)
{
	// Att!! Hi = \sum_j Pij/Pi log(Pij/Pi); (no negative sign !!)
	// because:
	// Pij = exp(-Bi *Lij)
	// log(Pij) = -Bi *Lij (note the sign !!!)
	Zi = -1.0;
	double Hi = .0;
	// +++ entropy (up to qNN)
	for (unsigned int j = 0; j < qNN; j++) {
		double Kij = std::exp(-Bi *Li[j]);
		Zi += Kij;
		Hi += Li[j] *Kij;
	}
	// Att!! Two problems may arise that cause Zi -> 0:
	// 1. outliers for which the starting value Bi = 1 is already too large
	// 2. for large Bi and largest values of Li[j], Kij might vanish and so Zi
	// 3. all qNN neighbours at the same distance makes Bi increase untill
	// they all fall out of the neighborhood making Zi = 0
	if (Zi > qNN *1e-12) {
		Hi *= Bi /Zi;
		Hi += std::log(Zi);
	}
	// In any case, Hi = 0 --> diff > 0 --> Bi is decreased
	else {
		Hi = 0;
	}
	return Hi;
}

// compute chunk of Betas
// [[Rcpp::export]]
Rcpp::NumericMatrix zBeta(int thread_rank, int threads, SEXP sexpX, bool is_distance, bool is_sparse, int ppx, double xppx)
{
	// input data
	sqDist* sqDistX = new sqDist(sexpX);
	// get chunk breaks
	std::vector<int> breaks (threads +1, 0);
	for (unsigned int b = 0; b < threads; b++) {
		breaks[b] = (int) b *(sqDistX ->nX +1.0) /threads;
	}
	breaks[threads] = sqDistX ->nX;
	// output data (chunk of Betas and theta-quantiles)
	int chunk_size = breaks[thread_rank] - breaks[thread_rank -1];
	Rcpp::NumericMatrix Beta(3, chunk_size);
	// parameters
	double logppx = std::log(ppx);
	double ppxtol = std::log(1.0 /.99999);
	double lambert_factor = -6.28318530718 /std::pow((xppx *ppx), 2.0);
	// quantile 3*ppx (L_max to compute Beta_i)
	int q3p = std::min((int) (3 *ppx +1), (int) (sqDistX ->nX -1));
	// quantile xppx*ppx (NN set boundary)
	int qNN = std::min((int) (xppx *ppx +1), (int) (sqDistX ->nX -1));
	//
	for (unsigned int i = breaks[thread_rank -1], l = 0; i < breaks[thread_rank]; i++, l++) {
		// squared distances
		std::vector<double> Li(sqDistX ->nX);
		if (is_sparse) {
			sqDistX ->row_spDist(i, Li.data());
		}
		else if (is_distance) {
			sqDistX ->row_d2Dist(i, Li.data());
		}
		else {
			sqDistX ->row_d1Dist(i, Li.data());
		}
		// nearest-neighbor quantile
		std::nth_element(Li.begin(), Li.begin() +qNN, Li.end());
		// compute Beta
		double B_lambertWm1 = -1.0 /(2 *Li[qNN]) *lambertWm1_CS(lambert_factor *Li[qNN]);
		double Bi = 1.0, minBeta = 0, maxBeta = DBL_MAX;
		double Zi;
		for (unsigned int iter = 0; iter < 100; iter++) {
			// compute log-perplexity for current Beta
			double Hi = entropy(Li, Bi, Zi, qNN);
			if ((Hi == 0) && (Bi > B_lambertWm1)) {
				Bi = B_lambertWm1;
				entropy(Li, Bi, Zi, qNN);
				break;
			}
			double diff = logppx -Hi;
			// double diff = logppx -Hi;
			// double Ppx = exp(-Bi *Li[ppx]);
			// printf("+++ %3u, %6.4e %6.4e %6.4e, %6.4e %6.4e %6.4e \n", i, Bi, Zi, Hi, B_lambertWm1, Ppx, Ppx /Zi);
			// check tolerance
			if (std::abs(diff) < ppxtol) break;
			// adjust Beta
			if (diff < 0){
				minBeta = Bi;
				if (maxBeta == DBL_MAX) Bi *= 2.0;
				else Bi = (Bi +maxBeta) /2.0;
			}
			else {
				maxBeta = Bi;
				Bi = (Bi +minBeta) /2.0;
			}
		}
		// printf("+++ %3u, %6.4e \n", i, Bi);
		// precision
		Beta(0, l) = Bi;
		// conditional affinity normalization factor (compute up to qNN)
		// std::nth_element(Li.begin(), Li.begin() +qNN, Li.end());
		// entropy(Li, Bi, Zi, qNN);
		// I drop the 1/n factor because affinities are renormalized
		// by thread (by factor zP) when attractive forces are computed
		Beta(1, l) = Zi *2.0;
		// nearest-neighbor quantile
		Beta(2, l) = Li[qNN];
	}
	// free memory
	delete sqDistX;
	//
	return Beta;
}

// Alternative method. LOOKS LIKE IT DOES NOT WORK !!
// entropy (up to qNN) (my version)
double entropy_lambertWm1(std::vector<double> Li, double Bi, double &Zi, double &Lx, int qNN)
{
	// Att!! Hi = \sum_j Pij/Pi log(Pij/Pi); (no negative sign !!)
	// because:
	// Pij = exp(-Bi *Lij)
	// log(Pij) = -Bi *Lij (note the sign !!!)
	Zi = -1.0;
	Lx = .0;
	double Hi = .0;
	// +++ entropy (up to qNN)
	for (unsigned int j = 0; j < qNN; j++) {
		double Kij = std::exp(-Bi *Li[j]);
		Zi += Kij;
		Hi += Li[j] *Kij;
		Lx += Li[j] *Kij;
	}
	// Att!! Two problems may arise that cause Zi -> 0:
	// 1. outliers for which the starting value Bi = 1 is already too large
	// 2. for large Bi and largest values of Li[j], Kij might vanish and so Zi
	// 3. all qNN neighbours at the same distance makes Bi increase untill
	// they all fall out of the neighborhood making Zi = 0
	Lx /= Zi;
	if (Zi > qNN *1e-12) {
		Hi *= Bi /Zi;
		Hi += std::log(Zi);
	}
	// In any case, Hi = 0 --> diff > 0 --> Bi is decreased
	else {
		Hi = 0;
	}
	return Hi;
}

// compute chunk of Betas
// [[Rcpp::export]]
Rcpp::NumericMatrix zBeta_lambertWm1(size_t thread_rank, size_t threads, SEXP sexpX, bool is_distance, bool is_sparse, size_t ppx, double xppx)
{
	// input data
	sqDist* sqDistX = new sqDist(sexpX);
	// get chunk breaks
	std::vector<int> breaks (threads +1, 0);
	for (size_t b = 0; b < threads; b++) breaks[b] = (int) b *(sqDistX ->nX +1.0) /threads;
	breaks[threads] = sqDistX ->nX;
	// output data (chunk of Betas and theta-quantiles)
	size_t chunk_size = breaks[thread_rank] - breaks[thread_rank -1];
	Rcpp::NumericMatrix Beta(3, chunk_size);
	// parameters
	double logppx = std::log(ppx);
	double ppxtol = std::log(1 /.99999);
	int qNN = std::min((int) (xppx *ppx +1), (int) (sqDistX ->nX -1));
	double lambert_factor = -6.28318530718 /std::pow((xppx *ppx), 2.0);
	//
	for (unsigned int i = breaks[thread_rank -1], l = 0; i < breaks[thread_rank]; i++, l++) {
		// squared distances
		std::vector<double> Li(sqDistX ->nX);
		if (is_sparse) {
			sqDistX ->row_spDist(i, Li.data());
		}
		else if (is_distance) {
			sqDistX ->row_d2Dist(i, Li.data());
		}
		else {
			sqDistX ->row_d1Dist(i, Li.data());
		}
		// nearest-neighbor quantile
		std::nth_element(Li.begin(), Li.begin() +qNN, Li.end());
		// starting (minimum) beta
		double Bi = -1.0 /(2 *Li[qNN]) *lambertWm1_CS(lambert_factor *Li[qNN]);
		// expected L
		double Lx;
		double Zi;
		// compute Beta
		for (size_t iter = 0; iter < 100; iter++) {
			// compute log-perplexity for current Beta
			// double diff = logppx -entropy_lambertWm1(Li, Bi, Zi, Lx, qNN);
			double Hi = entropy_lambertWm1(Li, Bi, Zi, Lx, qNN);
			double diff = logppx -Hi;
			double Ppx = exp(-Bi *Li[ppx]);
			printf("+++ %3u, %6.4e %6.4e %6.4e, %6.4e %6.4e %6.4e \n", i, Bi, Zi, Hi, Li[qNN], Lx, Ppx /Zi);
			// check tolerance
			if (std::abs(diff) < ppxtol) break;
			// adjust Beta
			Bi = -1.0 /(2 *Lx) *lambertWm1_CS(lambert_factor *Lx);
		}
		printf("+++ %3u, %6.4e \n", i, Bi);
		// precision
		Beta(0, l) = Bi;
		// conditional affinity normalization factor
		// I can drop the 1/n factor because affinities are renormalized
		// by thread (by factor zP) when attractive forces are computed
		// Beta(1, l) = Zi *2.0 *sqDistX ->nX;
		Beta(1, l) = Zi *2.0;
		// nearest-neighbor quantile
		Beta(2, l) = Li[qNN];
	}
	// free memory
	delete sqDistX;
	//
	return Beta;
}

/*** R

# Att!! After many checks with P1M, MNIST, S3D, Technosperm, Macosko15, ..
# I find out the settings below seem to work well:
# 1. Zi > qNN *1e-12	(this is the most definitive !!!)
# 2. logppx = log(ppx)
# 3. q3p = qNN = xppx *ppx (xppx = threads /layers)
# 4. Bi = 1
# 5. there is no point in center/scale the input data.

check1 <- function(D, ppx = 100, threads = 100, thread_rank = 1, is.distance = F, is.sparse = F, xppx = 3)
{
	Xbm <- bigmemory::as.big.matrix(t(D), type = 'double')
	# centerScale(Xbm@address, is.distance, is.sparse)
	zBeta(thread_rank, threads, Xbm@address, is.distance, is.sparse, ppx, xppx)
}

check2 <- function(D, ppx = 100, threads = 100, thread_rank = 1, is.distance = F, is.sparse = F, xppx = 3)
{
	Xbm <- bigmemory::as.big.matrix(t(D), type = 'double')
	zBeta_lambertWm1(thread_rank, threads, Xbm@address, is.distance, is.sparse, ppx, xppx)
}

scale <- function(D, is.distance, is.sparse)
{
	Xbm <- bigmemory::as.big.matrix(t(D), type = 'double')
	centerScale(Xbm@address, is.distance, is.sparse)
	return(t(Xbm[, ]))
}
*/
