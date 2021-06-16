/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <Rcpp.h>

#include "affmtx.h"
#include "qtree.h"
#include "tsne.h"

// ony for X2P dgemm_ version !!
// extern "C" {
//     #include <R_ext/BLAS.h>
// }

using namespace std;

// t-SNE constructor
TSNE::TSNE(unsigned int z, unsigned int w, unsigned int mY, double* eRange, int max_iter, double lRate, double theta, double alpha, double zP, double* exgg) : z(z), w(w), mY(mY), eRange(eRange), max_iter(max_iter), lRate(lRate), theta(theta), alpha(alpha), zP(zP), exgg(exgg)
{
	// set current value for learning rate
	std::vector<double> eta (mY, 0);
	for (unsigned int d = 0; d < mY; d++) eta[d] = lRate *4.0;
	this ->eta = eta;
}

// Perform t-SNE (only for 2D embedding !!)
void TSNE::run2D(double* P, int* W, double* Y)
{
	// +++ Allocate memory
	// . gradient forces
	double* atrF = (double*) calloc(z *mY, sizeof(double));
	double* repF = (double*) calloc(z *mY, sizeof(double));
	// . update of mapped positions
	double* uY = (double*) calloc(z *mY, sizeof(double));
	// +++ optimization
	for (int iter = 0; iter < max_iter; iter++) {
		double zQ = .0;
		if (theta < .33) {
			zQ = exact_Gradient(P, Y, atrF, repF);
		} else {
			zQ = apprx_Gradient(P, W, Y, atrF, repF);
		}
		// update (with momentum and learning rate)
		for (unsigned int i = 0, k = 0; i < z; i++){
			for (unsigned int d = 0; d < mY; d++, k++){
				// update embedding position
				uY[k] = alpha *uY[k] -eta[d] *(exgg[i] *atrF[k] /zP -repF[k] /zQ);
				Y[k] += uY[k];
				// reset attractive/repulsive forces
				atrF[k] = .0;
				repF[k] = .0;
			}
		}
	}
	// eRange is used to compute the embedding size (i.e. pooled from the master)
	// if not used to compute the learning-rate
	// there is no need to push it from the master (remove pushing in bdm_ptsne.R !!)
	for (unsigned int i = 0, k = 0; i < z; i++){
		for (unsigned int d = 0; d < mY; d++, k++){
			if (Y[k] < eRange[d *mY +0]) eRange[d *mY +0] = Y[k];
			else if (Y[k] > eRange[d *mY +1]) eRange[d *mY +1] = Y[k];
		}
	}
	// +++ Clean up memory
	free(atrF); atrF  = NULL;
	free(repF); repF = NULL;
	free(uY); uY = NULL;
	return;
}

// Compute gradient forces of the t-SNE cost function (EXACT)
double TSNE::exact_Gradient(double* P, double* Y, double* atrF, double* repF)
{
	double zQ = .0;
	double L [mY];
	for(unsigned int i = 0, ij = 0; i < z; i++) {
		for (unsigned int j = i +1; j < z; j++, ij++) {
			double Lij = 1.0;
			for(unsigned int d = 0; d < mY; d++) {
				L[d] = Y[i *mY +d] -Y[j *mY +d];
				Lij += L[d] *L[d];
			}
			double Qij = 1.0 /Lij;
			for(unsigned int d = 0; d < mY; d++) {
				if (P[ij] > 0) {
					double Aij = P[ij] *Qij *L[d];
					atrF[i *mY +d] += Aij;
					atrF[j *mY +d] -= Aij;
				}
				double Rij = Qij *Qij *L[d];
				repF[i *mY +d] += Rij;
				repF[j *mY +d] -= Rij;
			}
			zQ += Qij;
		}
	}
	zQ *= 2.0;
	return zQ;
}

// Compute gradient forces of the t-SNE cost function (APPRX)
double TSNE::apprx_Gradient(double* P, int* W, double* Y, double* atrF, double* repF)
{
	// compute attractive forces
	double L [mY];
	for(unsigned int k = 0; k < w; k++) {
		unsigned int i = W[k] /z;
		unsigned int j = (int) W[k] %z;
		double Lij = 1.0;
		for(unsigned int d = 0; d < mY; d++){
			L[d] = Y[i *mY +d] - Y[j *mY +d];
			Lij +=  L[d] *L[d];
		}
		unsigned int ij = ijIdx(z, i, j);
		for(unsigned int d = 0; d < mY; d++) {
			double Aij = P[ij] *L[d] /Lij;
			atrF[i *mY +d] += Aij;
			atrF[j *mY +d] -= Aij;
		}
	}
	// compute repulsive forces
	double zQ = .0;
	Quadtree* qtree = new Quadtree(Y, z, mY);
	qtree->repForces(Y, repF, theta, &zQ);
	delete qtree;
	//
	return zQ;
}

// Evaluate t-SNE cost function (exactly)
double TSNE::getCost(double* P, double* Y)
{
	// thread cost function:
	// Cost = -\sum_{i!=j} p_ij /zP * log (q_ij /zQ)
	// Q normalization factor
	double zQ = .0;
	// join cross entropy
	double jxH = .0;
	for (unsigned int i = 0, ij = 0; i < z; i++) {
		for (unsigned int j = i +1; j < z; j++, ij++){
			double Lij = 1.0;
			for(unsigned int d = 0; d < mY; d++) Lij += std::pow(Y[i *mY +d] -Y[j *mY +d], 2);
			if (P[ij] > 0) {
				jxH -= P[ij] *std::log(Lij);
			}
			zQ += 1.0 /Lij;
		}
	}
	// pseudo-normalized Cost
	jxH *= 2.0 /zP;
	jxH -= std::log(2.0 *zQ);
	double Cost = -jxH /std::log(z *(z -1));
	return Cost;
}
