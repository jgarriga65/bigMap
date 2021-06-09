/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

#ifndef TSNE_H
#define TSNE_H

class TSNE
{
public:
	size_t mY;							// embedding dimensions
	size_t z;							// thread size
	size_t w;							// approx. gradient size
	double* eRange;						// embedding range
	int max_iter;						// t-SNE parameters
	double lRate, theta, alpha, zP;
	std::vector<double> eta;
	double* exgg;
	// constructor
	TSNE(size_t z, size_t w, size_t mY, double* eRange, int max_iter, double lRate, double theta, double alpha, double zP, double* exgg);
	// run
	void run2D(double* P, int* W, double* Y);
	// compute Cost
	double getCost(double* P, double* Y);
	// update eRange
	void row_Gradient(double* P, int* W, double* Y, double* thread_Y, size_t i, size_t zi);
private:
	// compute gradient forces
	double exact_Gradient(double* P, double* Y, double* atrF, double* repF);
	double apprx_Gradient(double* P, int* W, double* Y, double* atrF, double* repF);
};

#endif
