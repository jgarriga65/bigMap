/*

+++ Structure divergence is computed by means of the Hellinger distance H(P,Q), a normalized symetric distance measure between probability distributions, given as,

	H(P,Q) = 1 /sqrt(2) srqt(\sum_i (sqrt(p_i) -sqrt(q_i))**2) = ... = sqrt(1 -\sum_i sqrt(p_i *q_i))
	(*1)

	herein, structure preservation can be defined as,

	SP = 1 - H(P, Q) = 1 - sqrt(1 -\sum_i sqrt(p_i *q_i))

+++ Local/Global stucture can be computed separately by considering the conditional distributions pairs (p_{j|i}, q_{j|i}) and (p_{i|j}, q_{i|j}) respectively. Then the overall local structure preservation is defined as the expected value of point-wise local structure preservation, that is,

	SP_i = 1 - H(P_{j|i}, Q_{j|i}) = 1 - sqrt(1 -\sum_j sqrt(p_j|i *q_j|i))

	SP = \sum_i P_{i} SP_{i} = \sum_i P_i * (1 - sqrt(1 -\sum_j sqrt(p_j|i *q_j|i)))
		= 1 - \sum_i P_i * sqrt(1 -\sum_j sqrt(p_j|i *q_j|i))

	and analogously for the global structure preservation.

+++ Alternatively, we can compute the expected value of the Hellinger distance,

	SD = \sum_{i} P_{i} H(P_{j|i}, Q_{j|i}) = \sum_{i} P_{i} sqrt(1 -\sum_j sqrt(p_j|i *q_j|i))

	and afterwards, compute the structure preservation as:

	SP = 1 - SD = 1 - \sum_{i} P_{i} sqrt(1 -\sum_j sqrt(p_j|i *q_j|i))

	what yields to the same expression as before.

+++ Att1!! Remember the following relation,

	(sqrt(p_i) - sqrt(q_i))**2 = (sqrt(p_i) + sqrt(q_i))**2 - 4 *sqrt(p_i*q_i)

	to understand the geometric interpretation of this measure

+++ Att2!! Equation (*1) is true AS LONG AS \sum_{i} p_i = 1 and \sum_{i} q_i = 1. As we can not be sure that this is so, it is better to relay on the original definition, that is,

 	H(P,Q) = 1/sqrt(2) sqrt(\sum_i (sqrt(p_i) -sqrt(q_i))**2)

+++ In any case, the marginals p_{i} are computed as,

	p_{i} = sum_{j} p_{ij} = sum_{j} (p_{j|i} +p_{i|j})/(2n) = 1/(2n) (1 + \sum_{j} p_{i|j})

	Att3!! Again, I'm assuming here that for all i \sum_{j} p_{j|i} = 1, what might not be true in practice.

*/

#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <math.h>
#include <float.h>
#include <vector>       // std::vector
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort, std::set_intersection

#include "mydist.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

static inline size_t tIdx(size_t n, size_t i, size_t j) {
	return  (n*i) - (i+1)*(i+2)/2 + j ;
}


// structure preservation: sp = 1 - H(P, Q)**2 (Hellinger distance)
// [[Rcpp::export]]
arma::Col<double> z_spQlty(SEXP sexpX, SEXP sexpB, SEXP sexpY, const arma::Col<int>& I, bool is_distance)
{
	// input data
	Rcpp::XPtr<BigMatrix> bmX(sexpX);
	MatrixAccessor<double> X(*bmX);
	// input Betas
	Rcpp::XPtr<BigMatrix> bmB(sexpB);
	MatrixAccessor<double> B(*bmB);
	// output data
	Rcpp::XPtr<BigMatrix> bmY(sexpY);
	MatrixAccessor<double> Y(*bmY);
	//
	size_t n = I.size();
	size_t m = bmX -> ncol();
	size_t d = bmY -> ncol();

	// P marginal distributions. Init 1 comes from \sum_j(p_{j|i}) (read below)
	std::vector<double> pM(n, 1.0);
	// P (NOT-normalized) conditional distribution
	std::vector<double> pC(n *n, .0);
	// P conditional distribution normalization factors (by rows and by columns)
	std::vector<double> pRow(n, .0);
	std::vector<double> pCol(n, .0);

	// Q (NOT-normalized) conditional distribution
	std::vector<double> qC(n *n, .0);
	// Q conditional distribution normalization factors (by rows and by columns)
	std::vector<double> qRow(n, .0);
	std::vector<double> qCol(n, .0);

	double Lij ;
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (j != i){
				if (is_distance) {
					Lij = X[I[j]][I[i]] * X[I[j]][I[i]];
				} else {
					Lij = .0;
					for(size_t v = 0; v < m; v++) Lij += std::pow(X[v][I[i]] - X[v][I[j]], 2);
				}
				pC[i *n +j] = std::exp(-B[0][I[i]] *Lij);
				pRow[i] += pC[i *n +j];
				pCol[j] += pC[i *n +j];
				Lij = .0;
				for(size_t v = 0; v < d; v++) Lij += std::pow(Y[v][I[j]] - Y[v][I[i]], 2);
				qC[i *n +j] = 1.0 / (1.0 + Lij);
				qRow[i] += qC[i *n +j];
				qCol[j] += qC[i *n +j];
			}
		}
		if (pRow[i] > 0) {				// read Att3!! above
			for (size_t j = 0; j < n; j++) {
				if (j != i) pM[j] += pC[i *n +j] /pRow[i];
			}
		} else pM[i] -= 1.0;			// compensate initialization value
	}

	//
	double lSD = .0;					// local structure divergence
	double gSD = .0;					// global structure divergence
	double sumM = .0;
	for (size_t i = 0; i < n; i++) {
		double lSDi = .0;
		if (pRow[i] > 0) {
			for (size_t j = 0; j < n; j++){
				if (j != i) {
					size_t ji = i *n +j;
					double pCji = pC[ji] /pRow[i];
					double qCji = qC[ji] /qRow[i];
					// lSDi += std::sqrt(pCji *qCji);	// read Att2!!! above
					lSDi += std::pow(std::sqrt(pCji) - std::sqrt(qCji), 2);
				}
			}
		} else lSDi = 2;
		double gSDi = .0;
		if (pCol[i] > 0) {
			for (size_t j = 0; j < n; j++){
				if (j != i) {
					size_t ij = j *n +i;
					// Att.! normalization is for column i here
					double pCij = pC[ij] /pCol[i];
					double qCij = qC[ij] /qCol[i];
					// gSDi += std::sqrt(pCij *qCij);	// read Att2!!! above
					gSDi += std::pow(std::sqrt(pCij) - std::sqrt(qCij), 2);
				}
			}
		} else gSDi = 2;
		// complete normalization of marginal P
		pM[i] /= (2 *n);
		// divergence expectation
		// lSD += pM[i] *std::sqrt(1 -lSDi);	// read Att2!!! above
		// gSD += pM[i] *std::sqrt(1 -gSDi);
		lSD += pM[i] *std::sqrt(lSDi) /std::sqrt(2);
		gSD += pM[i] *std::sqrt(gSDi) /std::sqrt(2);
		// this is to re-normalize the marginals in case we have some pRow[i] = 0
		sumM += pM[i];
	}

	// quality index
	arma::Col<double> Qlty(4);
	Qlty[0] = 1.0 - lSD /sumM;
	Qlty[1] = .0;
	Qlty[2] = 1.0 - gSD /sumM;
	Qlty[3] = .0;

	return Qlty;
}


// Rank-based Index of structure preservation
// [[Rcpp::export]]
arma::Col<double> z_rbQlty(SEXP sexpX, SEXP sexpY, const arma::Col<int>& I, bool is_distance)
{
	// input data
	Rcpp::XPtr<BigMatrix> bmX(sexpX);
	MatrixAccessor<double> X(*bmX);
	// output data
	Rcpp::XPtr<BigMatrix> bmY(sexpY);
	MatrixAccessor<double> Y(*bmY);
	//
	size_t n = I.size();
	size_t m = bmX -> ncol();
	size_t d = bmY -> ncol();
	// n = 10;

	// HD/LD local affinities
	std::vector<double> lP(n *n, .0);
	std::vector<double> lQ(n *n, .0);
	// HD/LD global affinities
	std::vector<double> gP(n *n, .0);
	std::vector<double> gQ(n *n, .0);

	std::vector<double> Li(n, 0);
	std::vector<int> rank(n);
	int k = 0;

	// +++ compute HD/LD rank-based local affinities
	for (size_t i = 0; i < n; i++) {

		// compute HD local distances
		for (size_t j = 0; j < n; j++) {
			if (is_distance) {
				Li[j] = X[I[j]][I[i]];
			} else {
				Li[j] = .0;
				for(size_t v = 0; v < m; v++) Li[j] += std::abs(X[v][I[i]] - X[v][I[j]]);
			}
		}
		// sort HD local distances
		k = 0;
		std::iota(rank.begin(), rank.end(), k ++);
 		std::stable_sort(rank.begin(), rank.end(), [&](int a, int b) {
			return Li[a] < Li[b]; });
		// compute HD local rank-based affinities
		lP[i *n +rank[0]] = .0;
		for (size_t j = 1; j < n; j++) lP[i *n +rank[j]] = 2.0 *(n -j) /(n *(n -1));
		// use gP to save Li[j] (by columns) for later computation of HD global affinities
		for (size_t j = 0; j < n; j++) gP[rank[j] *n +i] = j;

		// compute LD local distances
		for (size_t j = 0; j < n; j++) {
			Li[j] = .0;
			for(size_t v = 0; v < d; v++) Li[j] += std::abs(Y[v][I[i]] - Y[v][I[j]]);
		}
		// sort LD local distances
		k = 0;
		std::iota(rank.begin(), rank.end(), k ++);
 		std::stable_sort(rank.begin(), rank.end(), [&](int a, int b) {
			return Li[a] < Li[b]; });
		// compute LD local rank-based affinities
		lQ[i *n +rank[0]] = .0;
		for (size_t j = 1; j < n; j++) lQ[i *n +rank[j]] = 2.0 *(n -j) /(n * (n -1));
		// use gQ to save Li[j] (by columns) for later computation of LD global affinities
		for (size_t j = 0; j < n; j++) gQ[rank[j] *n +i] = j;

	}

	// +++ compue HD/LD rank-based global affinities
	for (size_t i = 0; i < n; i++) {
		// sort HD global distances
		k = 0;
		std::iota(rank.begin(), rank.end(), k ++);
 		std::stable_sort(rank.begin(), rank.end(), [&](int a, int b) {
			return gP[i *n +a] < gP[i *n +b]; });
		// compute HD global rank-based affinities
		gP[i *n +rank[0]] = .0;
		for (size_t j = 1; j < n; j++) gP[i *n +rank[j]] = 2.0 *(n -j) /(n *(n -1));
		// sort LD global distances
		k = 0;
		std::iota(rank.begin(), rank.end(), k ++);
 		std::stable_sort(rank.begin(), rank.end(), [&](int a, int b) {
			return gQ[i *n +a] < gQ[i *n +b]; });
		// if (thread_rank == 0) {
		// 	printf("+++ %2zu,", i);
		// 	for (size_t j = 0; j < n; j++) printf(" %2d    ", rank[j]);
		// 	printf("\n");
		// }
		// compute HD global rank-based affinities
		gQ[i *n +rank[0]] = .0;
		for (size_t j = 1; j < n; j++) gQ[i *n +rank[j]] = 2.0 *(n -j) /(n *(n -1));
		// if (thread_rank == 0) {
		// 	printf("+++ %2zu,", i);
		// 	for (size_t j = 0; j < n; j++) printf(" %4.4f", gQ[i *n +rank[j]]);
		// 	printf("\n");
		// }
	}

	double lRec = .0;						// local structure recall
	double lPrc = .0;						// local structure precision
	double gRec = .0;						// global structure recall
	double gPrc = .0;						// global structure precision
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++){
			if (j != i) {
				size_t ji = i *n +j;
				lRec += lP[ji] *std::log(lP[ji] /lQ[ji]);
				lPrc += lQ[ji] *std::log(lQ[ji] /lP[ji]);
				gRec += gP[ji] *std::log(gP[ji] /gQ[ji]);
				gPrc += gQ[ji] *std::log(gQ[ji] /gP[ji]);
			}
		}
	}

	arma::Col<double> Qlty(4);
	Qlty[0] = 1 -gRec /(n *std::log(n -1));
	Qlty[1] = 1 -gPrc /(n *std::log(n -1));
	Qlty[2] = 1 -lRec /(n *std::log(n -1));
	Qlty[3] = 1 -lPrc /(n *std::log(n -1));

	return Qlty;
}


// k-ary neighbourhood preservation
// [[Rcpp::export]]
arma::Mat<int> z_kNP(int thread_rank, int threads, SEXP sexpX, SEXP sexpY, bool is_distance, bool is_sparse, const arma::Col<int>& K, double sampling)
{
	// input data
	sqDist* sqDistX = new sqDist(sexpX);
	// output data
	Rcpp::XPtr<BigMatrix> bmY(sexpY);
	MatrixAccessor<double> Y(*bmY);
	int d = bmY -> ncol();
	int n = bmY -> nrow();
	// size of the set of k-ary-neighbourhoods
	int K_size = K.size();
	// get chunk breaks
	size_t chunkBeg = (size_t) (thread_rank -1) *(n +1.0) / threads;
	size_t chunkEnd = (size_t) (thread_rank == threads) ? n : thread_rank *(n +1.0) / threads;
	// sample data chunk (do not check all data-points for large datasets)
	std::vector<int> chunkIdx;
	size_t chunkSize = chunkEnd -chunkBeg;
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
		// squared distances
		std::vector<double> Li(sqDistX ->nX);
		if (is_sparse) {
			sqDistX ->row_spDist(*i, Li.data());
		}
		else if (is_distance) {
			sqDistX ->row_d2Dist(*i, Li.data());
		}
		else {
			sqDistX ->row_d1Dist(*i, Li.data());
		}
		// sort HD local distances
		std::vector<int> hd_rank(n);
		k = 0;
		std::iota(hd_rank.begin(), hd_rank.end(), k++);
 		std::stable_sort(hd_rank.begin(), hd_rank.end(), [&](int a, int b) {
			return Li[a] < Li[b]; });
		// compute LD local distances
		for (int j = 0; j < n; j++) {
			Li[j] = .0;
			for(int v = 0; v < d; v++) Li[j] += std::pow(Y[v][*i] - Y[v][j], 2);
		}
		// sort LD local distances
		std::vector<int> ld_rank(n);
		k = 0;
		std::iota(ld_rank.begin(), ld_rank.end(), k++);
 		std::stable_sort(ld_rank.begin(), ld_rank.end(), [&](int a, int b) {
			return Li[a] < Li[b]; });
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
	for (uword k = 0; k < zSampleSize; k++) {
		// datapoint index
		uword i = zSample[k] *sqDistX ->nX;
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
