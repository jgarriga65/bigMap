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

	unsigned int mY;					// embedding dimensions
	unsigned int z;						// thread size
	unsigned int nnSize;				// nearest neighbors size

	int max_iter;
	double theta, eta, alpha, gain;
	double zP;
	
	// constructor
	TSNE(unsigned int z, unsigned int nnSize, unsigned int mY, int max_iter, double theta, double lRate, double alpha, double gain, double zP);

	// run
	void run2D(double* P, unsigned int* W, double* Y);
	// compute Cost
	double Cost(double* P, unsigned int* W, double* Y);

private:
	// compute gradient forces
	void Gradient(double* P, unsigned int* W, double* Y, double* atrF, double* repF, double& zQ);
};

#endif
