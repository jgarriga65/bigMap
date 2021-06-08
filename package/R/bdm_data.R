# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Process raw input data

# Att!!
# Max. object size in R (< 3.0.0) is 2^31, thus maxRows * MaxCols:
#  10000 * 214748
# 100000 * 21478
#     1M * 2147
#    10M * 214

bdm.data <- function(raw.data, whiten = 4, input.dim = NULL, is.distance = F, is.sparse = F, quiet = TRUE)
{
	if (!quiet) cat('+++ processing data \n')
	if (is.null(input.dim)) {
		if (is.distance) {
			input.dim <- 1
		} else {
			input.dim <- ncol(raw.data)
		}
	}
	if (class(raw.data) == 'character') {
		return(list(inp.data = raw.data, whiten = 0, input.dim = input.dim))
	}
	else {
		if (whiten == 0){
			if (class(raw.data) != 'big.matrix') {
				return(list(inp.data = as.big.matrix(raw.data[, 1:input.dim], type = 'double'), whiten = 0, input.dim = input.dim))
			} else {
				return(list(inp.data = raw.data, whiten = 0, input.dim = input.dim))
			}
		}
		else {
			return(data.get(raw.data, whiten, input.dim))
		}
	}
}


# -----------------------------------------------------------------------------
# +++ Preprocessing of input-data
# -----------------------------------------------------------------------------

data.get <- function(raw.data, whiten, input.dim)
{
	# filter out all irrelevant features
	# (makes sense in datasets like the mnist.optical.digits where some features might have zeros for all observations)
	# X <- X[, which(apply(X, 2, sum) != 0)]
	# input.dim <- min(input.dim, ncol(X))

	if (whiten == 1)	# centering
	{
		X <- scale(raw.data[ , ], center = T, scale = F)
	}
	else if (whiten == 2)	# centering & scaling
	{
		X <- scale(raw.data[ , ], center = T, scale = T)
		if (any(is.na(X))) {
			return(message('+++ Error: scaling return NaNs !!!'))
		}
	}
	else if (whiten == 3 || whiten == 4)	# PCA/whitening
	{
		# TODO: consider using bigalgebra::bigPCA()
		X <- t(scale(raw.data[ , ], center = T, scale = F))
		# covariance matrix
		# Att!! ncol(X) stands for nrow(t(X)), the original X
		V <- X %*% t(X) / (ncol(X) -1)
		# singular value decomposition
		s <- svd(V, nu = input.dim, nv = 0)
		# PCA
		K <- t(s$u)
		# whitening
		if (whiten == 4) K <- diag(1/sqrt(s$d[1:input.dim])) %*% K
		# take first input.dim dimensions
		X <- t(K %*% X)
	}

	return(list(inp.data = as.big.matrix(X, type = 'double'), whiten = whiten, input.dim = input.dim))
}
