# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' Initialize MPI parallel computing environment.
#'
#' @param threads The number of parallel threads (in principle only limited by hardware resources, \code{i.e.} number of cores and available memory)
#'
#' @return cl A cluster instance (as created by the snow::makeCluster() function).

bdm.mpi.start <- function(threads)
{
	cl <- snow::makeCluster()
	# assert we have 1 core per thread plus 1 core (thread 0) running the master process
	if (Rmpi::mpi.comm.size(comm = 0) != (threads +1)) {
		message('+++ Error, for MPI clusters set cores = threads +1 \n')
		stopCluster(cl)
		cl <- NULL
	}
	else
	{
		cat('+++ starting ', threads, ' threads \n', sep = '')
		# cluster structure
		cl.ranks <- unlist(clusterEvalQ(cl, thread.rank <- Rmpi::mpi.comm.rank(comm = 0)))
		cl.nodes <- unlist(clusterEvalQ(cl, thread.node <- Rmpi::mpi.get.processor.name()))
		cl.hldrs <- sapply(unique(cl.nodes), function(node) min(cl.ranks[which(cl.nodes == node)]))
		# assign holders to workers
		clusterExport(cl, c('cl.hldrs'), envir = environment())
		nulL <- clusterEvalQ(cl, thread.hldr <- cl.hldrs[Rmpi::mpi.get.processor.name()])
	}
	return(cl)
}

#' Stops MPI parallel computing environment.
#'
#' @param cl A cluster instance (as created by the bdm.mpi.start() function).

bdm.mpi.stop <- function(cl)
{
	if (!is.null(cl)) stopCluster(cl)
}

bdm_.mpi.export <- function(cl, X, is.distance = F, is.sparse = F, normalize = T)
{
	Xdata.exp(cl, X, is.distance, is.sparse, normalize)
}

# Starts parallel computing environment.
#
# @return A cluster instance (as created by the snow::makeCluster() function).

cluster.start <- function(threads, mpi.cl = NULL, verbose = T)
{
	cl <- NULL
	if (!is.null(mpi.cl))
	{
		cl <- mpi.cl
	}
	else
	{
		if (verbose) cat('+++ starting ', threads, ' threads \n', sep = '')
		cl <- makeCluster(threads +1, type = 'PSOCK')
		# cluster structure
		# Att!!! bdm_glbl.R declares global variables: thread.rank, ...
		clusterApply(cl, seq_along(cl), function(i) thread.rank <<- i -1)
	}
	if (!is.null(cl))
	{
		clusterExport(cl, c('threads'), envir = environment())
		# load workers environment
		clusterEvalQ(cl, library(bigmemory))
		clusterEvalQ(cl, library(bigMap))
	} else if (verbose) {
		cat('+++ Error starting cluster !! \n')
	}
	# if (substr(bdm.local(), 1, 7) == 'xxx.xxx'){
	# 	if (verbose) cat('+++ WARNING: bdm.local() not set !! \n')
	# }
	return(cl)
}

# Close parallel computing environment.
#
# @param cl A cluster instance (as created by the snow::makeCluster() function).
#
# @return None

cluster.stop <- function(cl)
{
	if (!is.null(cl)) {
		if (attr(cl[[1]], 'class') == 'SOCKnode') stopCluster(cl)
	}
}
