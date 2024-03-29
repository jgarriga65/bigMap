# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' Clustering statistics box-plot.
#'
#' @param data A matrix of data to be plotted (either the input data matrix or any covariate matrix as well).
#'
#' @param bdm  A \var{bdm} data mapping instance.
#'
#' @param byVars Default value is \var{byVars = FALSE}, i.e. box-plots are grouped by cluster. If \var{byVars = TRUE} box-plots are grouped by input feature.
#'
#' @param merged Default value is \var{merged = TRUE}, i.e. the box-plot depicts the clusters after merging (if merging has been performed). If \var{merged = FALSE} or no merging has been performed the box-plot depicts the top-level clustering.
#'
#' @param clusters A vector with a subset of cluster ids. Default value is \var{clusters=NULL}, i.e. plot all clusters, (with a maximum of 25).
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @return None.
#'
#' @details If the number of clusters is large, only the first 25 clusters will be plotted. Note that the WTT algorithm numbers the clusters based on density value at the peak cell of the cluster. Thus, the numbering of the clusters is highly correlated with their relevance in terms of partial density.
#'
#' @examples
#'
#' bdm.example()
#' bdm.boxp(ex$data, ex$map)
#' bdm.boxp(ex$data, ex$map, byVars = TRUE)

bdm.boxp <- function(data, bdm, byVars = F, merged = T, clusters = NULL, layer = 1)
{
	# get data
	if (is.null(data))
		return(message('+++ Error: no data given !'))
	if (is.null(bdm$wtt)) {
		return(message('+++ Error: up-stream step WTT not found'))
	}

	bdm$data <- as.big.matrix(data, type = 'double')

	L <- bdm.labels(bdm, merged = merged, layer = layer)
	if (is.null(clusters)) {
		if (!is.null(bdm$merge))
			K <- sort(unique(bdm$merge$C))
		else
			K <- sort(unique(bdm$wtt[[layer]]$C))
	}
	else
		K <- clusters

	if (length(K) > 25) {
		message('+++ Showing only 25 first clusters')
		K <- K[1:25]
	}

	parbdm.set(mar = c(1.25, 1.25, 1.25, 0.75))

	if (byVars)
	{
		layout(layout.get(ncol(bdm$data)))
		nulL <- lapply(seq(ncol(bdm$data)), function(j) {
			X <- sapply(K, function(k) bdm$data[which(L == k), j])
			boxplot(X, main = colnames(bdm$data)[j], names = paste('C', K, sep = ''))
		})
		nullPlts <- max(layout.get(ncol(bdm$data))) -ncol(bdm$data) -1
	}
	else
	{
		layout(layout.get(length(K)))
		nulL <- lapply(K, function(k) {
			boxplot(bdm$data[which(L ==k), ], main = paste('cluster', k))
		})
		nullPlts <- max(layout.get(length(K))) -length(K) -1
	}

	if (nullPlts > 1) nulL <- sapply(seq(nullPlts-1), function(p) plot.null())

	parbdm.def()

}
