# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#' Find optimal number of clusters based on signal-to-noise-ratio.
#'
#' Performs a recursive merging of clusters based on minimum loss of signal-to-noise-ratio (S2NR). The S2NR is the explained/unexplained variance ratio measured in the high dimensional space based on the given low dimensional clustering. Merging is applied recursively until reaching a configuration of only 2 clusters and the S2NR is measured at each step.
#'
#' @param data Input data (a matrix, a big.matrix or a .csv file name).
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param plot.optk Default value is \var{plot.optk = TRUE}, i.e. the function plots the S2NR measure versus the number of clusters.
#'
#' @param ret.optk Default value is \var{ret.optk = FALSE}. For large datasets this computation can take a while, so it can be saved. If TRUE, the S2NR computation is saved and returned along with a \var{bdm} data mapping instance.
#'
#' @param info Default value is \var{info = TRUE}. If TRUE, all merging steps are shown.
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @return None if \var{ret.optk = FALSE}, else, a \var{bdm} data mapping instance.
#'
#' @details The underlying heuristic is that neigbouring clusters in the embedding correspond to close clusters in the high dimensional space, \var{i.e.} this merging heuristic is based on the spatial distribution of clusters. For each cluster (child cluster) we choose the neighboring cluster with steepest gradient along their common border (father cluster). Thus, we get a set of pairs of clusters (child/father) to be potentially merged. Given this set of candidates, the merging is performed recursively choosing, at each step, the pair of child/father clusters that results in a minimum loss of S2NR.
#' Typically some clusters dominate over all of their neighboring clusters. These clusters have no \var{father}. Thus, once all posible mergings have been performed we reach a \var{blocked} state where only the dominant clusters remain. This situation identifies a hierarchy level in the clustering. When this situation is reached, the algorithm starts a new merging round, identifying the child/father relations at that level of the hierarchy. The process stops when only two clusters remain.
#' Usually, the clustering hierarchy is clearly depicted by singular points in the S2NR function. This is a hint that the low dimensional clustering configuration is an image of a hierarchycal configuration in the high dimensional space. See \var{bdm.optk.plot()}.
#'
#' @examples
#'
#' # --- load mapped dataset
#' bdm.example()
#' # --- compute the optimal number of clusters and plot the S2NR
#' bdm.optk.s2nr(ex$data, ex$map, plot.optk = TRUE, ret.optk = FALSE)

bdm.optk.s2nr <- function(data, bdm, info = T, plot.optk = T, ret.optk = F, layer = 1)
{
	if (is.null(data))
		return(message('+++ Error: must define raw data !'))
	if (is.null(bdm$wtt[[layer]]$C)) {
		return(message('+++ Error: up-stream step bdm.wtt(layer = ', layer, ') not found !', sep = ''))
	}
	else {
		bdm$optk <- s2nr.optk(data, bdm, verbose = info, layer = layer)
		if (plot.optk) bdm.optk.plot(bdm)
		if (ret.optk) return(bdm)
	}
}


#' Plots the signal-to-noise-ratio as a function of the number of clusters.
#'
#' The function \var{bdm.optk.sn2r()} computes the S2NR that results from recursively merging clusters and, by deafult, makes a plot of these values. For large datasets this computation can take a while, so we can save this result by setting \var{optk.ret = TRUE}. If this result is saved, we can plot it again at any time using this funcion.
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @return None.
#'
#' @examples
#'
#' bdm.example()
#' m <- bdm.optk.s2nr(ex$data, ex$map, ret.optk = TRUE)
#' bdm.optk.plot(m)

bdm.optk.plot <- function(bdm)
{
	if (is.null(bdm$optk))
		return(message('+++ Error: up-stream step bdm.optk() is not computed ! \n'))

	optk <- bdm$optk
	s <- length(optk$H) +1
	# gradient function
	dffH <- c(diff(optk$H), 0)

	parbdm.set(oma = c(1,1,1,1), mar = c(3,3,2,0.5))
	layout(matrix(seq(2), c(2, 1)))

	plot(optk$H, xlab = '#clusters', ylab = 'S2NR', xaxt = 'n', type = 'l', col = '4')
	axis(1, at = seq(0, s, 10), labels = seq(s, 0, -10) +1, cex.axis = 0.6)
	grid()

	if (!is.null(optk$lvls))
	{
		L <- s -optk$lvls +1
		title('S2NR', adj = 0, cex = 0.7)
		abline(v = L, lwd = 2.0, col = "#BBBBBB")
		if (length(optk$lvls) < 4) {
			axis(3, at = L, labels = optk$lvls, cex.axis = 0.6, tick = FALSE)
		} else {
			sel <- seq(1, length(optk$lvls), by = 2)
			axis(3, at = L[sel], labels=optk$lvls[sel], cex.axis=0.6, tick=FALSE)
			sel <- seq(2, length(optk$lvls), by=2)
			axis(3, line=0.5, at=L[sel], labels=optk$lvls[sel], cex.axis=0.6, tick=FALSE)
		}
	}

	# if (!is.null(optk$tRate))
	# {
	# 	title('S2NR by transition.rate optimization', adj=0, cex=0.7)
	# 	points(optk$tRate, type='l', col='6')
	# 	legend("bottomleft", legend=c('S2NR', 'tRate'), col=c(4,6), cex=0.8, lwd=3, text.font=1, bty='n')
	# }

	plot(dffH, xlab = '#clusters', ylab = 'd(S2NR)', xaxt = 'n', type = 'l', col = '2')
	title('S2NR loss gradient', adj = 0, cex = 0.7)
	axis(1, at = seq(0, s, 10), labels = seq(s, 0, -10) +1, cex.axis = 0.6)
	grid()

	L <- s -optk$loss +1
	abline(v = L, lwd=2.0, col="#BBBBBB")
	if (length(optk$loss) < 4) {
		axis(3, at = L, labels = optk$loss, cex.axis = 0.6, tick = FALSE)
	} else{
		sel <- seq(1, length(optk$loss), by = 2)
		axis(3, at=L[sel], labels = optk$loss[sel], cex.axis = 0.6, tick = FALSE)
		sel <- seq(2, length(optk$loss), by = 2)
		axis(3, line = 0.5, at = L[sel], labels = optk$loss[sel], cex.axis = 0.6, tick = FALSE)
	}

	parbdm.def()
}

#' Merging of clusters based on signal-to-noise-ratio.
#'
#' Performs a recursive merging of clusters based on minimum loss of signal-to-noise-ratio (S2NR) until reaching the desired number of clusters. The S2NR is the explained/unexplained variance ratio measured in the high dimensional space based on the given low dimensional clustering.
#'
#' @param data Input data (a matrix, a big.matrix or a .csv file name).
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param k The number of desired clusters. The clustering will be recursively merged until reaching this number of clusters. Default value is \var{k = 10}. By setting \var{k < 0} we specify the number of clusters we want to merge (as opposed to the number of final clusters).
#'
#' @param plot.merge Default value is \var{plot.merge = TRUE}, i.e. the merged clustering is plotted.
#'
#' @param ret.merge Default value is \var{ret.merge = FALSE}. If \var{ret.merge = TRUE}, the function returns a \var{bdm} data mapping instance with the merged clustering.
#'
#' @param info Default value is \var{info = FALSE}. If TRUE, all merging steps are shown.
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @param ... If \var{plot.merge} is TRUE, you can pass plotting parameters to \var{bdm.wtt.plot()}.
#'
#' @return None if \var{ret.merge = FALSE}, else, a \var{bdm} data mapping instance.
#'
#' @details See details in \var{bdm.optk.s2nr()}.
#'
#' @examples
#'
#' # --- load mapped dataset
#' bdm.example()
#' \dontrun{
#' # --- merge and plot (no save)
#' bdm.merge.s2nr(ex$data, ex$map, k = 12)
#' # --- merge and return merging (do no plot)
#' m <- bdm.merge.s2nr(ex$data, ex$map, k = 12, ret.merge = T, plot.merge = F)
#' # --- plot merging afterwards (use bdm.wtt.plot() function)
#' bdm.wtt.plot(m)
#' }

bdm.merge.s2nr <- function(data, bdm, k = 10, plot.merge = T, ret.merge = F, info = T, layer = 1, ...){

	if (is.null(data))
		return(message('+++ Error: must define raw data !'))
	if (k >= bdm$wtt[[layer]]$s)
		return(message('+++ Error: k must be lower than the current number of clusters !'))
	bdm$merge <- s2nr.merge(data, bdm, k = k, verbose = info, layer = layer)
	if (length(bdm$merge$steps) == 0)
		return(message('+++ Merging error !!'))
	# show clustering after merge
	if (plot.merge) bdm.wtt.plot(bdm, layer = layer, ...)
	# return new hdd
	if (ret.merge) return(bdm)
}


# -----------------------------------------------------------------------------
# +++ Merge auxiliary functions
# -----------------------------------------------------------------------------

merge.labels <- function(Y, C, G){
	D2c <- grid_D2cell(Y, G) +1
	#g.size <- G[1, 1]
	lbls <- C[(D2c[, 2] - 1) *G[1, 1] +D2c[, 1]]
	return(lbls)
}
