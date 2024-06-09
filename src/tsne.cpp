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
#include <fstream>

#include "affmtx.h"
#include "qtree.h"
#include "tsne.h"

using namespace std;

// DBL_EPSILON = 1.0e-9 (DBL_MIN = 1.0e-37)

// SHOOT OFF CONDITION CONTROL
// . minL1
// static const double minAtrForce = 1.0e-6 /(1.0 +2.0e-12);
// .minL2 (this one makes sense but does not work either)
// static const double minAtrForce = std::sqrt(DBL_EPSILON) /(1.0 +2.0 *DBL_EPSILON);
// .minL3
// static const double minAtrForce = 2.0 *std::sqrt(DBL_EPSILON) /(1.0 +2.0 *DBL_EPSILON);
// .minL interval
// static const double lowL = 1.0 /std::sqrt(6.0);
// static const double uppL = 1.0 /std::sqrt(2.0);

// t-SNE constructor
TSNE::TSNE(unsigned int z, unsigned int nnSize, unsigned int mY, int max_iter, double theta, double lRate, double alpha, double gain, double zP) : z(z), nnSize(nnSize), mY(mY), max_iter(max_iter), theta(theta), eta(lRate), alpha(alpha), gain(gain), zP(zP)
{}

// Perform t-SNE (only for 2D embedding !!)
void TSNE::run2D(double* P, unsigned int* W, double* Y, int dbgRow, unsigned int thread_rank, unsigned int epoch)
{
	// . gradient forces
	double* atrF = (double*) calloc(z *mY, sizeof(double));
	double* repF = (double*) calloc(z *mY, sizeof(double));

	// . position updates
	double* updY = (double*) calloc(z *mY, sizeof(double));
	// . position gains
	double* etaG = (double*) calloc(z *mY, sizeof(double));

	// reset step size at each new epoch !!
	for (unsigned int k = 0; k < z *mY; k++) etaG[k] = 1.0;
	// double maxGain = 4.0;

	// . update clipping value
	// double uClip = eta /(2.0 *std::log(z *nnSize));

	// debugging option
	// char fName [12];
	// int lName = std::sprintf(fName, "epoch%02u_%03u", epoch, thread_rank);
    // std::ofstream dbgFile;
    // dbgFile.open(fName, std::ios::app);

	// +++ optimization
	double zQ = .0;
	for (int iter = 0; iter < max_iter; iter++) {

		TSNE::getForces(P, W, Y, atrF, repF, zQ);

		for (unsigned int i = 0, k = 0; i < z; i++) {
			for (unsigned int d = 0; d < mY; d++, k++) {

				// update (with momentum and learning rate)
				double dY =  atrF[k] /zP -repF[k] /zQ;
				if (gain > 0) {
					if (signbit(dY) != signbit(updY[k])) etaG[k] += (gain *.1);
					else etaG[k] *= (1.0 -gain /10.0);
					// if ((maxGain > 0) && (etaG[k] > maxGain)) etaG[k] = maxGain;
				}
				// update embedding position
				updY[k] = alpha *updY[k] -eta *etaG[k] *dY;
				Y[k] += updY[k];

				// debugging
				// if (i == dbgRow) {
				// 	dbgFile << epoch << ",";
				// 	dbgFile << iter << ",";
				// 	dbgFile << d << ",";
				// 	dbgFile << zP << ",";
				// 	dbgFile << atrF[k] /zP << ",";
				// 	dbgFile << zQ << ",";
				// 	dbgFile << repF[k] /zQ << ",";
				// 	dbgFile << dY << ",";
				// 	dbgFile << etaG[k] << ",";
				// 	dbgFile << updY[k] << ",";
				// 	dbgFile << Y[k] << "\n";
				// }

				// clip position update to the size of the embedding (shoot off control)
				// if (std::abs(Y[k]) > uClip){
				// 	Y[k] = std::signbit(Y[k]) ? -uClip : uClip;
				// 	updY[k] = 0.0;
				// 	etaG[k] = 1.0;
				// }

				// reset attractive/repulsive forces
				atrF[k] = .0;
				repF[k] = .0;
			}
		}
	}

	// debugging option
	// dbgFile.flush();
	// dbgFile.close();

	// +++ Clean up memory
	free(atrF); atrF  = NULL;
	free(repF); repF = NULL;
	free(updY); updY = NULL;
	free(etaG); etaG = NULL;
}

// Compute gradient forces of the t-SNE cost function
void TSNE::getForces(double* P, unsigned int* W, double* Y, double* atrF, double* repF, double& zQ)
{
	double L [mY];
	// attractive forces
	for(unsigned int i = 0; i < z; i++) {
		for (unsigned int ni = 0; ni < nnSize; ni++) {
			if (P[i *nnSize +ni] > 0) {
				unsigned int j = W[i *nnSize +ni];
				double Lij = 1.0;
				for(unsigned int d = 0; d < mY; d++) {
					L[d] = Y[i *mY +d] -Y[j *mY +d];
					Lij += (L[d] *L[d]);
				}
				double Qij = 1.0 /Lij;
				for(unsigned int d = 0; d < mY; d++) {
					atrF[i *mY +d] += P[i *nnSize +ni] *Qij *L[d];
				}
			}
		}
	}
	// repulsive forces
	zQ = .0;
	if (theta < 0.33) {
		// exact repulsive forces
		for(unsigned int i = 0; i < z; i++) {
			for (unsigned int j = i +1; j < z; j++) {
				double Lij = 1.0;
				for(unsigned int d = 0; d < mY; d++) {
					L[d] = Y[i *mY +d] -Y[j *mY +d];
					Lij += (L[d] *L[d]);
				}
				double Qij = 1.0 /Lij;
				for(unsigned int d = 0; d < mY; d++) {
					double Qd = Qij *L[d];
					repF[i *mY +d] += Qij *Qd;
					repF[j *mY +d] -= Qij *Qd;
				}
				zQ += Qij;
			}
		}
		zQ *= 2.0;
	}
	else {
		// apprx. repulsive forces
		Quadtree* qtree = new Quadtree(Y, z, mY, theta);
		qtree->repForces(Y, repF, &zQ);
		delete qtree;
	}
}

// Evaluate t-SNE cost function (exactly)
double TSNE::Cost(double* P, unsigned int* W, double* Y)
{
	// thread cost function:
	// Cost = -\sum_{i!=j} p_ij /zP * log (q_ij /zQ)
	// join cross entropy
	double jxH = .0;
	for (unsigned int i = 0; i < z; i++) {
		for (unsigned int ni = 0; ni < nnSize; ni++) {
			if (P[i *nnSize +ni] > 0) {
				unsigned int j = W[i *nnSize +ni];
				double Lij = 1.0;
				for(unsigned int d = 0; d < mY; d++) Lij += std::pow(Y[i *mY +d] -Y[j *mY +d], 2);
				jxH -= P[i *nnSize +ni] *std::log(Lij);
			}
		}
	}
	//
	double zQ = .0;
	for(unsigned int i = 0; i < z; i++) {
		for (unsigned int j = i +1; j < z; j++) {
			double Lij = 1.0;
			for(unsigned int d = 0; d < mY; d++) Lij += std::pow(Y[i *mY +d] -Y[j *mY +d], 2.0);
			zQ += 1.0 /Lij;
		}
	}
	// pseudo-normalized Cost
	jxH /= zP;
	jxH -= std::log(2.0 *zQ);
	double Cost = -jxH /std::log(z *(z -1));
	return Cost;
}
