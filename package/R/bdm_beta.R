# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.

# -----------------------------------------------------------------------------
# +++ Parallelized computation of perplexity-based local betas (by chunks)
# -----------------------------------------------------------------------------

# Compute perplexity-based local betas. This function assumes that \var{Xbm} (input-data big.matrix) is attached to the workers of \var{cl}. Also \var{is.distance} must have been exported to the workers. The up-stream step \code{bdm.data()} works out both conditions.

# Att!!
# beta_i = 1/(2 *sigma**2) and sigma_i = 1/sqrt(2 *beta_i)
# thus, the BIVARIATE normal density function:
# pdf(dij) = 1/(2 *pi *sigma_i^2) exp(-1/2 * dij^2/sigma_i^2)
# is writen in terms of betak as:
# pdf(dij) = 1/pi *beta_i exp(- beta_i * dij^2)

beta.get <- function(cl, ppx, xppx = 3.0)
{
	cat('+++ computing Betas, perplexity ', ppx, ' \n', sep='')
	t <- system.time({
		# export parameters
		clusterExport(cl, c('ppx', 'xppx'), envir=environment())
		# get perplexity-based local betas by chunks
		B <- matrix(unlist(clusterCall(cl, thread.beta)), ncol = 3, byrow = T)
		print(summary(B[, 1]))
	})
	print(t)
	return(list(ppx = ppx, xppx = xppx, B = B, t = t))
}

# -----------------------------------------------------------------------------
# +++ worker function
# -----------------------------------------------------------------------------

thread.beta <- function()
{
	if (thread.rank != 0) {
		zBeta(thread.rank, threads, Xbm@address, is.distance, is.sparse, ppx, xppx)
	}
}
