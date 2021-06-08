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

#include "mydist.h"
#include "lambertW.h"

// using namespace arma;
using namespace Rcpp;
using namespace std;


// entropy (up to 32.0/Bi or up to nnQ)
double entropy(std::vector<double> Li, double Bi, double &Zi, int nnQ)
{
	// Att!! Hi = \sum_j Pij/Pi log(Pij/Pi);
	// where:
	// Pij = exp(-Bi *Lij)
	// log(Pij) = -Bi *Lij (note the sign !!!)
	Zi = -1.0;
	double Hi = .0;
	// // entropy (up to 32.0/Bi)
	// double maxL = 32.0 /Bi;
	// for (size_t j = 0; j < Li.size(); j++) {
	// 	if (Li[j] <= maxL) {
	// 		double Kij = std::exp(-Bi *Li[j]);
	// 		Zi += Kij;
	// 		Hi += Li[j] *Kij;
	// 	}
	// }
	// entropy (up to nnQ)
	for (size_t j = 0; j < nnQ; j++) {
		double Kij = std::exp(-Bi *Li[j]);
		Zi += Kij;
		Hi += Li[j] *Kij;
	}
	// Att!! Two problems may arise that cause Zi = 0:
	// 1. datapoints for which the starting value Bi = 1 is too large
	// (i.e, distances to all neighbors are greater than 32.0/Bi)
	// 2. datapoints with many nearest neighbours at the same distance for which,
	// increasing Bi suddenly places them all out of the neighborhood
	// In any case, Hi = 0 and Bi is decreased in the next step.
	if (Zi > 0) {
		Hi *= Bi /Zi;
		Hi += std::log(Zi);
	}
	return Hi;
}

// compute chunk of Betas
// [[Rcpp::export]]
Rcpp::NumericMatrix zBeta(size_t thread_rank, size_t threads, SEXP sexpX, bool is_distance, bool is_sparse, size_t ppx, double xppx)
{
	// input data
	sqDist* sqDistX = new sqDist(sexpX);
	// get chunk breaks
	std::vector<int> breaks (threads +1, 0);
	for (size_t b = 0; b < threads; b++) breaks[b] = (int) b *(sqDistX ->nX +1.0) /threads;
	breaks[threads] = sqDistX ->nX;
	// output data (chunk of Betas and theta-quantiles)
	size_t chunk_size = breaks[thread_rank] - breaks[thread_rank -1];
	Rcpp::NumericMatrix Beta(4, chunk_size);
	// parameters
	double logppx = std::log(ppx);
	double ppxtol = std::log(1 /.99999);
	int nnQ = std::min((int) (xppx *ppx +1), (int) (sqDistX ->nX -1));
	double stdppx = std::pow((std::log(std::sqrt(6.28318530718)) -logppx), 2.0);
	double sq3ppx = std::pow(3.0 *ppx, 2.0);
	//
	for (size_t i = breaks[thread_rank -1]; i < breaks[thread_rank]; i++) {
		// chunk-index
		size_t l = i -breaks[thread_rank -1];
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
		std::nth_element(Li.begin(), Li.begin() +nnQ, Li.end());
		double nnQBeta = stdppx /Li[nnQ];
		// compute Beta
		double Bi = 1.0, minBeta = 0, maxBeta = DBL_MAX;
		double Zi;
		for (size_t iter = 0; iter < 100; iter++) {
			if (Bi > nnQBeta) {
				Bi = -1.0 /(2.0 *Li[nnQ]) *lambertWm1_CS(-(6.28318530718 *Li[nnQ]) /sq3ppx);
				// Bi = nnQBeta;
				entropy(Li, Bi, Zi, nnQ);
				break;
			}
			// compute log-perplexity for current Beta
			double diff = logppx -entropy(Li, Bi, Zi, nnQ);
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
		// precision
		Beta(0, l) = Bi;
		// conditional affinity normalization factor
		// I can drop the 1/n factor because affinities are renormalized
		// by thread (by factor zP) when attractive forces are computed
		// Beta(1, l) = Zi *2.0 *sqDistX ->nX;
		Beta(1, l) = Zi *2.0;
		// nearest-neighbor quantile
		Beta(2, l) = Li[nnQ];
		// exaggeration factor
		Beta(3, l) = 1.0;
		// Beta(3, l) = 12.0;
		// Beta(3, l) = 1.566555 /std::atan(4.5 /Bi);
		// Beta(3, l) = 159.551211035 /(std::sqrt(Bi) *3.141592 *(1 +4.5 /Bi));
		// Beta(3, l) = 159.335837 /(std:sqrt(Bi) *(1 + 4.5 /Bi) *(1.57079633 +std::atan(3 /std::sqrt(2 *Bi))));
	}
	// free memory
	delete sqDistX;
	//
	return Beta;
}
