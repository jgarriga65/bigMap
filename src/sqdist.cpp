/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "bigmemory/BigMatrix.h"
// [[Rcpp::depends(BH, bigmemory)]]

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "sqdist.h"

using namespace std;

double spDist(unsigned int m, double* Xa, double* Xb)
{
	double L = .0;
	unsigned int a = 0;
	unsigned int b = 0;
	while ((a < m && Xa[a] > 0) || (b < m && Xb[b] > 0))
	{
		if (b == m || Xb[b] == 0) {
			L += std::pow(Xa[a +1], 2);
			a += 2;
		}
		else if (a == m || Xa[a] == 0) {
			L += std::pow(Xb[b +1], 2);
			b += 2;
		}
		else if (Xa[a] < Xb[b]) {
			L += std::pow(Xa[a +1], 2);
			a += 2;
		}
		else if (Xb[b] < Xa[a]) {
			L += std::pow(Xb[b +1], 2);
			b += 2;
		}
		else if (Xa[a] == Xb[b]) {
			L += std::pow(Xa[a +1] -Xb[b +1], 2);
			a += 2;
			b += 2;
		}
	}
	return L;
}

sqDist::sqDist(SEXP sexpX)
{
	BigMatrix *bigX = reinterpret_cast <BigMatrix*> (R_ExternalPtrAddr(sexpX));
	index_type offset = bigX ->nrow() *bigX ->col_offset();
	this ->X = reinterpret_cast <double*> (bigX ->matrix()) +offset;
	// Att.!!! X transposed
	this ->nX = bigX ->ncol();
	this ->mX = bigX ->nrow();
	//
	// this ->minL = DBL_MAX;

	// Att.!!! És MOLT PERILLOS FER AIXÒ !!!!!
	// Aquesta matriu és compartida i per tan cada thread està executant això simultàniament.
	// Això només ho ha d'executar el thread.rank = 0 (implementat en la funció Xdata.exp())
	// scaleCenter(this ->X, this ->nX, this ->mX);

	bigX = NULL;
}

void sqDist::row_spDist(unsigned int i, double* Li)
{
	double Xi[this ->mX];
	double Xj[this ->mX];
	// this ->minL = DBL_MAX;
	for (unsigned int v = 0; v < this ->mX; v++) Xi[v] = this ->X[i *this ->mX +v];
	for (unsigned int j = 0; j < this ->nX; j++) {
		for (unsigned int v = 0; v < this ->mX; v++) Xj[v] = this ->X[j *this ->mX +v];
		Li[j] = spDist(mX, Xi, Xj);
	}
}

void sqDist::row_d2Dist(unsigned int i, double* Li)
{
	// this ->minL = DBL_MAX;
	for (unsigned int j = 0; j < this ->mX; j++) {
		Li[j] = std::pow(this ->X[i *this ->mX +j], 2);
	}
}

void sqDist::row_d1Dist(unsigned int i, double* Li)
{
	double Xi[this ->mX];
	// this ->minL = DBL_MAX;
	for (unsigned int v = 0; v < this ->mX; v++) Xi[v] = this ->X[i *this ->mX +v];
	for (unsigned int j = 0; j < this ->nX; j++) {
		Li[j] = .0;
		for(unsigned int v = 0; v < this ->mX; v++) Li[j] += std::pow(Xi[v] -this ->X[j *this ->mX +v], 2);
	}
}
