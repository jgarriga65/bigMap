# +++ The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.

# +++

#' HD/LD correlation.
#'
#' Pair-wise distance correlation between HD and LD neighborhoods.
#'
#' @param data Input data (a matrix, a big.matrix or a .csv file name).
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param zSampleSize Number of data points to check by thread. Default value is \var{zSampleSize = 1000}.
#'
#' @param mpi.cl An MPI (inter-node parallelization) cluster as returned by \var{bdm.mpi.start()}. Default value is \var{mpi.cl = NULL}, i.e. a 'SOCK' (intra-node parallelization) cluster is generated.
#'
#' @param threads The number of parallel threads (according to data size and hardware resources, i.e. number of cores and available memory). Default value is \var{threads = 4}.
#'
#' @return A \var{bdm} data mapping instance.
#'
#'
#' @examples
#'
#' # --- load example dataset
#' \dontrun{
#' bdm.example()
#' m <- bdm.hlCorr(ex$data, ex$map, threads = 4)
#' }
#'

bdm.hlCorr <- function(data, bdm, zSampleSize = 1000, mpi.cl = NULL, threads = 4)
{
	cl <- cluster.start(threads, mpi.cl)
	# . input data (if using mpi.cl it might have been already exported)
	if (is.null(mpi.cl) || !is.null(data)) {
		cat('+++ exporting input data \n')
		Xdata.exp(cl, data, bdm$is.distance, bdm$is.sparse, bdm$normalize)
	}
	cat('+++ exporting output data \n')
	Ydata.exp(cl, t(bdm$ptsne$Y[, 1:2]))
	cat('+++ computing hl-Correlation \n')
	clusterExport(cl, c('zSampleSize'), envir=environment())
	t <- system.time({
		bdm$hlC <- unlist(clusterCall(cl, thread.hlCorr))
	})
	print(summary(bdm$hlC))
	print(t)
	bdm$t$hlC <- t
	cluster.stop(cl)
	return(bdm)
}

thread.hlCorr <- function()
{
	if (thread.rank > 0)
		z_hlCorr(Xbm@address, Ybm@address, zSampleSize, is.distance, is.sparse)
}
