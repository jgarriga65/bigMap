/* lambertW.cpp

 Copyright (C) 2015, Avraham Adler
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

References:

Corless, R. M.; Gonnet, G. H.; Hare, D. E.; Jeffrey, D. J. & Knuth, D. E.
 "On the Lambert W function", Advances in Computational Mathematics,
 Springer, 1996, 5, 329-359

Fritsch, F. N.; Shafer, R. E. & Crowley, W. P.
 "Solution of the transcendental equation (we^w = x)",
 Communications of the ACM, Association for Computing Machinery (ACM),
 1973, 16, 123-124
*/

#ifndef LAMW_H
#define LAMW_H

const double EPS = 2.2204460492503131e-16;
const double M_1_E = 1.0 / M_E;

  /* Fritsch Iteration
  * W_{n+1} = W_n * (1 + e_n)
  * z_n = ln(x / W_n) - W_n
  * q_n = 2 * (1 + W_n) * (1 + W_n + 2 / 3 * z_n)
  * e_n = z_n / (1 + W_n) * (q_n - z_n) / (q_n - 2 * z_n)
  */

double FritschIter(double x, double w_guess){
	double w = w_guess;
	int MaxEval = 5;
	bool CONVERGED = false;
	double k = 2.0 / 3.0;
	int i = 0;
	do {
		double z = std::log(x / w) - w;
		double w1 = w + 1.0;
		double q = 2.0 * w1 * (w1 + k * z);
		double qmz = q - z;
		double e = z / w1 * qmz / (qmz - z);
		CONVERGED = std::abs(e) <= EPS;
		w *= (1.0 + e);
		++i;
	} while (!CONVERGED && i < MaxEval);
	return(w);
}


double lambertW0_CS(double x) {
	double result;
	double w;
	if (x == R_PosInf) {
		result = R_PosInf;
	} else if (x < -M_1_E) {
		result = R_NaN;
	} else if (std::abs(x + M_1_E) <= EPS) {
		result = -1.0;
	} else if (x <= M_E - 0.5) {
		if (std::abs(x) <= 1e-16) {
			/* This close to 0 the W_0 branch is best estimated by its Taylor/Pade
			expansion whose first term is the value x and remaining terms are below
			machine double precision. See
			https://math.stackexchange.com/questions/1700919/how-to-derive-the-lambert-w-function-series-expansion
			*/
			result = x;
		} else {
			if (std::abs(x) <= 7e-3) {
				/* Use equation (5) in Fritsch */
				w = ((1.33333333333333333 * x + 1.0) * x) /
				  ((0.83333333333333333 * x + 2.33333333333333333) * x + 1.0);
			} else {
				/* Use expansion in Corliss 4.22 to create (3, 2) Pade approximant
				Numerator:-10189 / 303840 * p^3 + 40529 / 303840 * p^2 + 489 / 844 * p-1
				Denominator: -14009 / 303840 * p^2 + 355 / 844 * p + 1
				Converted to digits to reduce needed operations
				*/
				double p = std::sqrt(2.0 * (M_E * x + 1.0));
				double Numer = ((-0.03353409689310163 * p + 0.1333892838335966) * p +
				                0.5793838862559242) * p - 1.0;
				double Denom = (-0.04610650342285413 * p + 0.4206161137440758) * p + 1.0;
				w = Numer / Denom;
			}
			result = FritschIter(x, w);
		}
	} else {
		/* Use first five terms of Corliss et al. 4.19 */
		w = std::log(x);
		double L_2 = std::log(w);
		double L_3 = L_2 / w;
		double L_3_sq = L_3 * L_3;
		w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
		w + L_3_sq * L_3 / 3.0;
		result = FritschIter(x, w);
	}
	return(result);
}

double lambertWm1_CS(double x){
	double result;
	double w;
	if (x == 0.0) {
		result = R_NegInf;
	} else if (x < -M_1_E || x > 0.0) {
		result = R_NaN;
	} else if (std::abs(x + M_1_E) <= EPS) {
		result = -1.0;
	} else {
		/* Use first five terms of Corliss et al. 4.19 */
		w = std::log(-x);
		double L_2 = std::log(-w);
		double L_3 = L_2 /w;
		double L_3_sq = L_3 *L_3;
		w += -L_2 +L_3 +0.5 *L_3_sq -L_3 /w +L_3 /(w *w) -1.5 *L_3_sq /w +L_3_sq *L_3 /3.0;
		result = FritschIter(x, w);
	}
	return(result);
}


#endif
