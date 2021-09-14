/*
*/

#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <math.h>
#include <float.h>
#include <vector>       // std::vector
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort, std::set_intersection

#include "sqdist.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

// k-ary neighbourhood preservation
// [[Rcpp::export]]
arma::Mat<int> z_kNP(int thread_rank, int threads, SEXP sexpX, SEXP sexpY, bool is_distance, bool is_sparse, const arma::Col<int>& K, double sampling)
{
	// input data
	sqDist* sqDistX = new sqDist(sexpX);
	// output data
	sqDist* sqDistY = new sqDist(sexpY);
	unsigned int n = sqDistX ->nX;
	// size of the set of k-ary-neighbourhoods
	int K_size = K.size();
	// get chunk breaks
	unsigned int chunkBeg = (unsigned int) (thread_rank -1) *(n +1.0) / threads;
	unsigned int chunkEnd = (unsigned int) (thread_rank == threads) ? n : thread_rank *(n +1.0) / threads;
	// sample data chunk (do not check all data-points for large datasets)
	std::vector<int> chunkIdx;
	unsigned int chunkSize = chunkEnd -chunkBeg;
	arma::Col<double> chunkRnd = arma::randu(chunkSize);
	for (int q = 0; q < chunkSize; q++) {
		if (chunkRnd[q] < sampling) chunkIdx.push_back(chunkBeg +q);
	}
	// output data
	arma::Mat<int> chunk_Q(chunkIdx.size(), K_size);
	//
	int k = 0;
	int q = 0;
	for (vector<int>::iterator i = chunkIdx.begin(); i != chunkIdx.end(); i++) {
		// HD pairwise squared distances
		std::vector<double> LHi(n);
		if (is_sparse) {
			sqDistX ->row_spDist(*i, LHi.data());
		}
		else if (is_distance) {
			sqDistX ->row_d2Dist(*i, LHi.data());
		}
		else {
			sqDistX ->row_d1Dist(*i, LHi.data());
		}
		// sort HD distances
		std::vector<int> hd_rank(n);
		k = 0;
		std::iota(hd_rank.begin(), hd_rank.end(), k++);
 		std::stable_sort(hd_rank.begin(), hd_rank.end(), [&](int a, int b) {
			return LHi[a] < LHi[b]; });
		// LD pairwise squared distances
		std::vector<double> LDi(n);
		sqDistY ->row_d1Dist(*i, LDi.data());
		// sort LD distances
		std::vector<int> ld_rank(n);
		k = 0;
		std::iota(ld_rank.begin(), ld_rank.end(), k++);
 		std::stable_sort(ld_rank.begin(), ld_rank.end(), [&](int a, int b) {
			return LDi[a] < LDi[b]; });
		// k-ary neihgborhood preservation
		for (int s = 0; s < K_size; s++) {
			// sort HD K[s] neighborhood by neighbor-number
			std::vector<int> hd_neigh(K[s]);
			std::copy(hd_rank.begin(), hd_rank.begin() +K[s], hd_neigh.begin());
			std::sort(hd_neigh.begin(), hd_neigh.end());
			// sort LD K[s] neighborhood by neighbor-number
			std::vector<int> ld_neigh(K[s]);
			std::copy(ld_rank.begin(), ld_rank.begin() +K[s], ld_neigh.begin());
			std::sort(ld_neigh.begin(), ld_neigh.end());
			// compute matching
			std::vector<int> match(2 *K[s]);
			std::vector<int>::iterator it;
			it = std::set_intersection(hd_neigh.begin(), hd_neigh.end(), ld_neigh.begin(), ld_neigh.end(), match.begin());
			match.resize(it -match.begin());
			// save
			chunk_Q(q, s) = match.size();
		}
		q++;
	}
	// free memory
	delete sqDistX;
	delete sqDistY;
	//
	return chunk_Q;
}

// Pairwise distance correlation
// [[Rcpp::export]]
double z_hlCorr(SEXP sexpX, SEXP sexpY, int zSampleSize, bool is_distance, bool is_sparse)
{
	// input data
	sqDist* sqDistX = new sqDist(sexpX);
	// output data
	sqDist* sqDistY = new sqDist(sexpY);
	// sample indexes
	arma::Col<double> zSample = arma::randu(zSampleSize);
	// correlation
	arma::Col<double> hlCorr(zSampleSize);
	for (unsigned int k = 0; k < zSampleSize; k++) {
		// datapoint index
		unsigned int i = zSample[k] *sqDistX ->nX;
		// HD pairwise distances
		arma::Col<double> HDi(sqDistX ->nX);
		if (is_sparse) {
			sqDistX ->row_spDist(i, HDi.memptr());
		}
		else if (is_distance) {
			sqDistX ->row_d2Dist(i, HDi.memptr());
		}
		else {
			sqDistX ->row_d1Dist(i, HDi.memptr());
		}
		// LD pairwise distances
		arma::Col<double> LDi(sqDistY ->nX);
		sqDistY ->row_d1Dist(i, LDi.memptr());
		// correlation
		arma::mat Ci = arma::cor(HDi, LDi, 0);
		hlCorr[k] = Ci(0, 0);
	}
	// free memory
	delete sqDistX;
	delete sqDistY;
	//
	return arma::mean(hlCorr);
}
