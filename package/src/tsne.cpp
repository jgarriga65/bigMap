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

using namespace std;

// t-SNE constructor
TSNE::TSNE(unsigned int z, unsigned int w, unsigned int mY, int max_iter, double theta, double* lRate, double alpha, double zP, double exgg) : z(z), w(w), mY(mY), max_iter(max_iter), theta(theta), alpha(alpha), zP(zP), exgg(exgg)
{
	// learning rate
	std::vector<double> eta (mY);
	for (size_t d = 0; d < mY; d++) eta[d] = lRate[d] *4.0;
	this ->eta = eta;
	// DBL_EPSILON = 1.0e-9 (DBL_MIN = 1.0e-37)
	this ->minL = DBL_EPSILON /4.0;
}

// Perform t-SNE (only for 2D embedding !!)
void TSNE::run2D(double* P, unsigned int* W, double* Y)
{
	// . gradient forces
	double* atrF = (double*) calloc(z *mY, sizeof(double));
	double* repF = (double*) calloc(z *mY, sizeof(double));
	// . mapping position updates
	double* updY = (double*) calloc(z *mY, sizeof(double));
	// . learning rate gains
	double* gain = (double*) calloc(z *mY, sizeof(double));
	for (unsigned int i = 0; i < z; i++) gain[i] = 1.0;
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
				// gradient
				double dY = exgg *atrF[k] /zP -repF[k] /zQ;
				// learning rate gain
				if (signbit(dY) != signbit(updY[k])) {
					gain[k] += .2;
				}
				else {
					gain[k] *= .8;
					gain[k] += .01;
				}
				// update embedding position
				updY[k] = alpha *updY[k] -eta[d] *gain[k] *dY;
				Y[k] += updY[k];
				// reset attractive/repulsive forces
				atrF[k] = .0;
				repF[k] = .0;
			}
		}
	}
	// +++ Clean up memory
	free(atrF); atrF  = NULL;
	free(repF); repF = NULL;
	free(updY); updY = NULL;
	free(gain); gain = NULL;
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
				if (std::abs(L[d]) < minL) {
					L[d] = (Y[i *mY +d] >Y[j *mY +d]) ? minL : -minL;
				}
				Lij += L[d] *L[d];
			}
			double Qij = 1.0 /Lij;
			for(unsigned int d = 0; d < mY; d++) {
				double Qd = Qij *L[d];
				if (P[ij] > 0) {
					atrF[i *mY +d] += P[ij] *Qd;
					atrF[j *mY +d] -= P[ij] *Qd;
				}
				repF[i *mY +d] += Qij *Qd;
				repF[j *mY +d] -= Qij *Qd;
			}
			zQ += Qij;
		}
	}
	zQ *= 2.0;
	return zQ;
}

// Compute gradient forces of the t-SNE cost function (APPRX)
double TSNE::apprx_Gradient(double* P, unsigned int* W, double* Y, double* atrF, double* repF)
{
	// compute attractive forces
	double L [mY];

	if (w < z) {
		// this is slower only if w << z
		for(unsigned int k = 0; k < w; k++) {
			unsigned int i = W[k] /z;
			unsigned int j = (unsigned int) W[k] %z;
			double Lij = 1.0;
			for(unsigned int d = 0; d < mY; d++){
				L[d] = Y[i *mY +d] - Y[j *mY +d];
				if (std::abs(L[d]) < minL) {
					L[d] = (Y[i *mY +d] >Y[j *mY +d]) ? minL : -minL;
				}
				Lij +=  L[d] *L[d];
			}
			unsigned int ij = ijIdx(z, i, j);
			for(unsigned int d = 0; d < mY; d++) {
				atrF[i *mY +d] += P[ij] *L[d] /Lij;
				atrF[j *mY +d] -= P[ij] *L[d] /Lij;
			}
		}
	}
	else {
		// otherwise (high perplexities) this is faster
		for(unsigned int i = 0, ij = 0; i < z; i++) {
			for (unsigned int j = i +1; j < z; j++, ij++) {
				if (P[ij] > 0) {
					double Lij = 1.0;
					for(unsigned int d = 0; d < mY; d++) {
						L[d] = Y[i *mY +d] -Y[j *mY +d];
						if (std::abs(L[d]) < minL) {
							L[d] = (Y[i *mY +d] >Y[j *mY +d]) ? minL : -minL;
						}
						Lij += L[d] *L[d];
					}
					for(unsigned int d = 0; d < mY; d++) {
						double Qd = L[d] /Lij;
						atrF[i *mY +d] += P[ij] *Qd;
						atrF[j *mY +d] -= P[ij] *Qd;
					}
				}
			}
		}
	}

	// compute repulsive forces
	double zQ = .0;
	Quadtree* qtree = new Quadtree(Y, z, mY, theta);
	qtree->repForces(Y, repF, &zQ);
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
