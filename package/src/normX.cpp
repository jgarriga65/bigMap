#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "bigmemory/BigMatrix.h"
// [[Rcpp::depends(BH, bigmemory)]]

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
	size_t nX = bigX ->ncol();
	size_t mX = bigX ->nrow();
	// +++ Center
	if (!is_distance)
	{
		std::vector<double> xMean(mX, 0);
		// if (is_sparse)
		// {
		// 	std::vector<double> xSize(mX, 0);
		// 	for (size_t i = 0, ij = 1; i < nX; i++) {
		// 		for (size_t j = 1; j < mX; j++, j++, ij++, ij++) {
		// 			xMean[j] += X[ij];
		// 			if (X[ij] > 0) xSize[j] ++;
		// 		}
		// 	}
		// 	// for (size_t j = 1; j < mX; j++, j++) xMean[j] /= (double) nX;
		// 	for (size_t j = 1; j < mX; j++, j++) xMean[j] /= xSize[j];
		// 	// Subtract data mean
		// 	for (size_t i = 0, ij = 1; i < nX; i++) {
		// 		for (size_t j = 1; j < mX; j++, j++, ij++, ij++) X[ij] -= xMean[j];
		// 	}
		// }
		// else
		if (!is_sparse)
		{
			for (size_t i = 0, ij = 0; i < nX; i++) {
				for (size_t j = 0; j < mX; j++, ij++) xMean[j] += X[ij];
			}
			for (size_t j = 0; j < mX; j++) xMean[j] /= (double) nX;
			// Subtract data mean
			for (size_t i = 0, ij = 0; i < nX; i++) {
				for (size_t j = 0; j < mX; j++, ij++) X[ij] -= xMean[j];
			}
		}
	}
	// +++ Scale
	double Xmax = .0;
	if (is_sparse)
	{
		for (size_t ij = 1; ij < nX *mX; ij++, ij++) if (std::abs(X[ij]) > Xmax) Xmax = std::abs(X[ij]);
		for (size_t ij = 1; ij < nX *mX; ij++, ij++) X[ij] /= Xmax;
	}
	else
	{
		for (size_t ij = 0; ij < nX *mX; ij++) if (std::abs(X[ij]) > Xmax) Xmax = std::abs(X[ij]);
		for (size_t ij = 0; ij < nX *mX; ij++) X[ij] /= Xmax;
	}
}
