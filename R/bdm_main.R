# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' Example dataset
#'
#' Loads a mapping example.
#'
#' @return An example dataset named \var{ex}
#'
#' @details The object \var{ex} is a list with elements: \var{ex$data}, a matrix with raw data; \var{ex$labels}, a vector of datapoint labels; \var{ex$map}, a \var{bdm} data mapping instance. A \var{bdm} instance is the basic object of the mapping protocol, i.e. a list to which new elements are added at each step of the mapping protocol.
#'
#' This example is based on a small synthetic dataset with \var{n = 5000} observations drawn from a 4-variate Gaussian Mixture Model (GMM) with 16 Gaussian components with random means and variances.
#'
#' @examples
#'
#' # --- load example dataset
#' bdm.example()
#' str(ex)

bdm.example <- function()
{
	load(paste(paste(system.file('extdata', package = 'bigMap'), '/', sep = ''), 'exMap.RData', sep=''), envir=parent.frame(), verbose=T)
}


#' Initialization of a \var{bdm} data mapping instance.
#'
#' Computes the precision parameters for the given perplexity (i.e. the local bandwidths for the input affinity kernels) and returns them as a \var{bdm} data mapping instance. A \var{bdm} data mapping instance is the starting object of the mapping protocol, a list to which new elements are added at each step of the mapping protocol.
#'
#' @param data A \var{data.frame} or \var{matrix} with raw input-data. The dataset must not have duplicated rows.
#'
#' @param is.distance Default value is \var{is.distance = FALSE}. TRUE indicates that raw data is a distance matrix.
#'
#' @param is.sparse Default value is \var{is.sparse = FALSE}. TRUE indicates that the raw data is a sparse matrix.
#'
#' @param ppx The value of perplexity to compute similarities.
#'
#' @param mpi.cl An MPI (inter-node parallelization) cluster as returned by \var{bdm.mpi.start()}. By default \var{mpi.cl = NULL}, i.e. a 'SOCK' (intra-node parallelization) cluster is used.
#'
#' @param threads Number of parallel threads (according to data size and hardware resources, \var{i.e.} number of cores and available memory. Default value is \var{threads = 4}).
#'
#' @param dSet.labels If available, labels can be included as a separate vector of length equal to \var{nrow(data)}. Label values are factorized as \var{as.numeric(as.factor(labels))}.
#'
#' @return A \var{bdm} data mapping instance. A \var{bdm} instance is the starting object of the mapping protocol, a list to which new elements are added at each step of the mapping protocol.
#'
#' @examples
#'
#' # --- load example dataset
#' bdm.example()
#' \dontrun{
#' m <- bdm.init(ex$data, ppx = 250, labels = ex$labels)
#' }

bdm.init <- function(data, is.distance = F, is.sparse = F, ppx = 100, mpi.cl = NULL, threads = 4, dSet.labels = NULL)
{
	bdm <- list()
	# check data is given as fileName.csv
	if (class(data) == 'character') {
	 	if (!file.exists(data)) {
			message('+++ Data file not found !! \n')
			return(FALSE)
		}
		bdm$dataFile <- data
	}
	bdm$is.distance <- is.distance
	bdm$is.sparse <- is.sparse
	bdm$normalize <- normalize
	# compute betas
	if (threads > 0) {
		# start cluster
		cl <- cluster.start(threads, mpi.cl)
		# export data
		cat('+++ exporting data \n')
		t <- system.time({
			Xdata.exp(cl, data, is.distance, is.sparse, normalize = normalize)
		})
		print(t)
		dimX <- Xdata.dim(cl)
		bdm$nX <- dimX[1]
		bdm$mX <- dimX[2]
		# compute betas
		if (length(ppx) == 1) {
			bdm$ppx <- list(beta.get(cl, ppx, xppx = xppx))
		} else {
			bdm$ppx <- lapply(ppx, function(ppx_i) beta.get(cl, ppx_i, xppx = xppx))
		}
		# stop cluster
		if (is.null(mpi.cl)) cluster.stop(cl)
	}
	# attach labels
	if (!is.null(dSet.labels)) bdm$lbls <- as.numeric(as.factor(dSet.labels))
	return(bdm)
}

#' Parallelized t-SNE (ptSNE)
#'
#' Starts the parallelized t-SNE algorithm (pt-SNE). This is the first step of the mapping protocol.
#'
#' @param data Input data (a matrix, a big.matrix or a .csv file name).
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param theta Accuracy/speed trade-off factor, a value between 0.33 and 0.8. (Default value is \var{theta = 0.0}). If \var{theta < 0.33} the algorithm uses the exact computation of the gradient. The closer is this value to 1 the faster is the computation but the coarser is the approximation of the gradient.
#'
#' @param Y.init A \var{n *2 *layers} matrix with initial mapping positions. (By default \var{Y.init=NULL} will use random initial positions).
#'
#' @param mpi.cl MPI (inter-node parallelization) cluster as generated by \var{bdm.mpi.start()}. (By default \var{mpi.cl = NULL} a 'SOCK' (intra-node parallelization) cluster is generated).
#'
#' @param threads Number of parallel threads (according to data size and hardware resources, \var{i.e.} number of cores and available memory. Default value is \var{threads = 4}).
#'
#' @param layers Number of layers (\var{minimum} 2, \var{maximum} the number of threads). Default value is \var{layers = 2}.
#'
#' @param info Output information: 1 yields inter-round results, 0 disables intermediate results. Default value is \var{info = 0}.
#'
#' @return A \var{bdm} data mapping instance.
#'
#' @examples
#'
#' # --- load example dataset
#' bdm.example()
#' # --- perform ptSNE
#' \dontrun{
#' # --- run ptSNE
#' m <- bdm.ptsne(ex$data, ex$map, threads = 10, layers = 2)
#' # --- plot the Cost function
#' bdm.cost(m)
#' # --- plot ptSNE output
#' bdm.ptsne.plot(m, class.lbls = ex$labels)
#' }

bdm.ptsne <- function(data, bdm, theta = 0.5, Y.init = NULL, mpi.cl = NULL, threads = 4, layers = 2, info = 0)
{
	# +++ sanity check
	if (!is.null(Y.init) && (ncol(Y.init) != 2 *layers)) {
		return(message('+++ ncol(Y.init) does not match the number of layers !! \n'))
	}
	if (layers > threads) {
		cat('+++ WARNING: layers set to ', threads, ' !!\n', sep='')
		layers <- threads
	}
	if (theta > 0.0 && theta < 0.33) {
		cat('+++ WARNING: theta set to ', 0.0, ' !!\n', sep='')
		theta <- 0.0
	} else if (theta > 0.8) {
		cat('+++ WARNING: theta set to ', 0.8, ' !!\n', sep='')
		theta <- 0.8
	}
	# +++ start cluster of workers
	bdm$t <- list()
	cl <- cluster.start(threads, mpi.cl)
	if (is.null(cl)) return(bdm)
	# export data (if using mpi.cl it might have been already exported)
	if (is.null(mpi.cl) || !is.null(data)) {
		cat('+++ exporting data \n')
		bdm$t$dataExport <- system.time({
			Xdata.exp(cl, data, bdm$is.distance, bdm$is.sparse, normalize = bdm$normalize)
		})
		print(bdm$t$dataExport)
	}
	#
	bdm$ptsne <- list(threads = threads, layers = layers, theta = theta, gain = gain, momentum = momentum, qDecay = qDecay, Y = Y.init)
	#
	bdm <- ptsne.get(cl, bdm, info)
	if (length(bdm) == 1) bdm <- bdm[[1]]
	# stop cluster
	cluster.stop(cl)
	return(bdm)
}

#' Restart pt-SNE
#'
#' Restarts the ptSNE algorithm (runs more epochs).
#'
#' @param data Input data (a matrix, a big.matrix or a .csv file name).
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param epochs Number of epochs to run. Default value \var{epochs = NULL} runs \var{4 *log(n)} epochs.
#'
#' @param iters Number of iters per epoch. Default value \var{iters = NULL} runs \var{4 *log(thread_size)} iters/epoch.
#'
#' @param mpi.cl An MPI (inter-node parallelization) cluster as returned by \var{bdm.mpi.start()}. Default value is \var{mpi.cl = NULL}, i.e. a 'SOCK' (intra-node parallelization) cluster is automatically generated.
#'
#' @param threads Number of parallel threads (according to data size and hardware resources, i.e. number of cores and available memory). Default value is \var{threads = 4}.
#'
#' @param layers Number of layers (\var{minimum} 2, \var{maximum} the number of threads). Default value is \var{layers = 2}.
#'
#' @param info Output information: 1 yields inter-round results, 0 disables intermediate results. Default value is 0.
#'
#' @return A \var{bdm} data mapping instance.
#'
#' @examples
#'
#' # --- load example dataset
#' bdm.example()
#' \dontrun{
#' # --- restart ptSNE
#' m <- bdm.restart(ex$data, ex$map, epochs = 50)
#' }

bdm.restart <- function(data, bdm, epochs = NULL, iters = NULL, mpi.cl = NULL, threads = NULL, layers = NULL, info = 0) {
	# +++ start cluster of workers
	if (!is.null(threads)) bdm$ptsne$threads <- threads
	if (!is.null(layers)) bdm$ptsne$layers <- layers
	cl <- cluster.start(bdm$ptsne$threads, mpi.cl)
	if (is.null(cl)) return(bdm)
	# export data (if using mpi.cl it might have been already exported)
	if (is.null(mpi.cl) || !is.null(data)) {
		cat('+++ exporting data \n')
		Xdata.exp(cl, data, bdm$is.distance, bdm$is.sparse, normalize = bdm$normalize)
	}
	# +++ run ptsne
	if (!is.null(epochs)) bdm$epochs <- epochs
	if (!is.null(iters)) bdm$iters <- iters
	cost <- bdm$ptsne$cost
	size <- bdm$ptsne$size
	bdm <- ptsne.restart(cl, bdm, info)
	if (!is.null(cost)) {
		if (bdm$ptsne$threads > 1) {
		 	if (nrow(cost) == nrow(bdm$ptsne$cost)) {
				bdm$ptsne$cost <- cbind(cost, bdm$ptsne$cost)
			}
		} else {
			bdm$ptsne$cost <- c(cost, bdm$ptsne$cost)
		}
	}
	bdm$ptsne$size <- c(size, bdm$ptsne$size)
	# +++ stop cluster
	cluster.stop(cl)
	return(bdm)
}

#' Perplexity-adaptive kernel density estimation
#'
#' Starts the paKDE algorithm (second step of the mapping protocol).
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param ppx The value of perplexity to compute similarities in the low-dimensional embedding. Default value is \var{ppx = 100}.
#'
#' @param g The resolution of the density space grid (\eqn{g*g} cells). Default value is \var{g = 200}.
#'
#' @param g.exp A numeric factor to avoid border effects. The grid limits will be expanded so as to enclose the density of the kernel of the most extreme embedded datapoints up to \var{g.exp} times \eqn{\sigma}. Default value is \var{g.exp = 3}, \var{i.e.} the grid limits are expanded so as to enclose the 0.9986 of the probability mass of the most extreme kernels.
#'
#' @param mpi.cl An MPI (inter-node parallelization) cluster as returned by \var{bdm.mpi.start()}. Default value is \var{mpi.cl = NULL}, i.e. a 'SOCK' (intra-node parallelization) cluster is automatically generated.
#'
#' @param threads Number of parallel threads (according to data size and hardware resources, i.e. number of cores and available memory). Default value is \var{threads = 4}.
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @details When computing the \var{paKDE} the embedding area is discretized as a grid of size \var{g*g} cells. In order to avoid border effects, the limits of the grid are expanded by default so as to enclose at least the 0.9986 of the cumulative distribution function (\eqn{3 \sigma}) of the kernels of the most extreme mapped points in each direction.
#'
#' The presence of outliers in the embedding can lead to undesired expansion of the grid limits. We can overcome this using lower values of \var{g.exp}. By setting \var{g.exp = 0} the grid limits will be equal to the range of the embedding.
#'
#' The values \var{g.exp = c(1, 2, 3, 4, 5, 6)} enclose cdf values of \var{0.8413, 0.9772, 0.9986, 0.99996, 0.99999, 1.0} respectively.
#'
#' @return A copy of the input \var{bdm} instance with new element \var{bdm$pakde} (paKDE output). \var{bdm$pakde[[layer]]$layer = 'NC'} stands for not computed layers.
#'

#' @examples
#'
#' # --- load mapped dataset
#' bdm.example()
#' # --- run paKDE
#' \dontrun{
#' m <- bdm.pakde(ex$map, ppx = 200, g = 200, g.exp = 3, threads = 4)
#' # --- plot paKDE output
#' bdm.pakde.plot(m)
#' }

bdm.pakde <- function(bdm, ppx = 100, g = 200, g.exp = 3, mpi.cl = NULL, threads = 2, layer = 1)
{
	# start cluster of workers
	cl <- cluster.start(threads, mpi.cl)
	if (is.null(cl)) return(bdm)
	# initialize bdm$pakde
	if (is.null(bdm$pakde)) bdm$pakde <- list()
	# compute kde
	l <- c(1, 2) + (layer- 1) *2
	if (!is.null(bdm$ptsne$Y[ , l])) {
		cat('+++ paKDE for layer ', layer, '/', bdm$ptsne$layers, ' +++ \n', sep='')
		bdm$pakde[[layer]] <- pakde.get(cl, bdm$ptsne$Y[ , l], ppx, g, g.exp)
	}
	else cat('+++ Error: up-stream step bdm.ptsne(layer = ', layer, ') not found ! \n', sep = '')
	# stop cluster
	cluster.stop(cl)
	return(bdm)
}


#' Watertrack transform (WTT)
#'
#' Starts the WTT algorithm (third setp of the mapping protocol).
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @return A \var{bdm} data mapping instance.
#'
#' @details This function requires the up-stream step \var{bdm.pakde()}.
#'
#' @examples
#'
#' # --- load mapped dataset
#' bdm.example()
#' # --- perform WTT
#' m <- bdm.wtt(ex$map)
#' # --- plot WTT output
#' bdm.wtt.plot(m)

bdm.wtt <- function(bdm, layer = 1)
{
	# initialize bdm$wtt
	if (is.null(bdm$wtt)) bdm$wtt <- list()
	# compute WTT
	if (!is.null(bdm$pakde[[layer]]$z))
	{
		cat('\n')
		cat('+++ WTT for layer ', layer, '/', bdm$ptsne$layers, ' +++ \n', sep='')
		bdm$wtt[[layer]] <- wtt.get(bdm$pakde[[layer]])
	}
	else cat('+++ Error: up-stream step bdm.pakde(layer = ', layer, ') not found ! \n', sep = '')
	return(bdm)
}


#' Get data-point clustering labels.
#'
#' Given that clusters are computed at grid-cell level, this function returns the clustering label for each data-point.
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @param merged Default value is \var{merged = TRUE}. If \var{merged = TRUE} and the clustering has been merged, the labels are the ids of the clusters after merging. If \var{merged = FALSE} or the clustering has not been merged, the labels indicate the ids of to the top-level clustering.
#'
#' @return A vector of data-point clustering labels.
#'
#' @examples
#'
#' bdm.example()
#' m.labels <- bdm.labels(ex$map)

bdm.labels <- function(bdm, merged = T, layer = 1){
	# At.!!! there is an internal version of this function (for simplicity)
	# check merge.labels() in bdm_merge.R if any change is to be made here !!
	if (!is.null(bdm$wtt[[layer]]))
	{
		C <- bdm$wtt[[layer]]$C
		if (merged && !is.null(bdm$merge)) {
			C <- bdm$merge$C
		}
		l <- c(1, 2) + (layer -1) *2
		D2c <- grid_D2cell(bdm$ptsne$Y[ , l], bdm$wtt[[layer]]$grid) +1
		x.size <- bdm$wtt[[layer]]$grid[1, 1]
		lbls <- C[(D2c[, 2] - 1) *x.size +D2c[, 1]]
	}
	else {
		lbls <- rep(1, nrow(bdm$ptsne$Y))
	}
	return(lbls)
}
