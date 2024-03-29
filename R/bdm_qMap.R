# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' ptSNE quantile-maps
#'
#' Overlay of quantitative variables onto the embedding.
#'
#' @param data A \var{matrix/data.frame} of data to overlay onto the embedding. This function generates a multi-plot layout with one overlay for each column in \var{data}.
#'
#' @param bdm A \var{bdm} data mapping instance.
#'
#' @param labels A vector of class labels of length \var{nrow(bdm$data)} to overlay onto the embedding. Label values are factorized as \var{as.numeric(as.factor(labels))}. Default value is \var{labels = NULL}.
#'
#' @param subset A numeric vector with the indexes of a subset of data. The subset of data-points is heat-mapped and the rest of the data-points are shown in light grey. Default is \var{subset = NULL}, i.e. all data-points are heat-mapped.
#'
#' @param qMap.levels The number of levels of the quantile-map. Default value is \var{qMap.levels = 8}.
#'
#' @param qMap.cex The size of the data-points (as in \var{par()}).
#'
#' @param qMap.bg The background colour of the qMap plot. Default value is \var{ptsne.bg = #FFFFFF} (white).
#'
#' @param class.pltt A colour palette to indicate class labels. Default value is \var{qMap.pltt = NULL}, i.e. a default palette is used.
#'
#' @param qtitle A vector of strings with titles for the plots. Default value is \var{qtitle = NULL}.
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @return None.
#'
#' @details This is not a heat-map but a quantile-map plot. This function splits the range of each variable into as many quantiles as specified by \var{levels} so that the color gradient will hardly ever correspond to a constant numeric gradient. Thus, the mapping will show more evenly distributed colors though at the expense of possibly exaggerating artifacts. For variables with very extrem distributions, it will be impossible to find as many quantiles as desired and the distribution of colors will not be so homogeneous.
#' @examples
#'
#' bdm.example()
#' bdm.qMap(ex$data, ex$map)
#' # --- show only components (1, 2, 4, 8) of the GMM
#' bdm.qMap(ex$data, ex$map, subset = which(ex$labels %in% c(1, 2, 4, 8)))

bdm.qMap <- function(data, bdm, labels = NULL, subset = NULL, qMap.levels = 8, qMap.cex = 0.3, qMap.bg = '#FFFFFF', class.pltt = NULL, qtitle =NULL, layer = 1)
{
	# get data
	if (is.null(data))
		return(message('+++ Error: no data given !'))

	# check 1 single plot
	if (class(data) == 'numeric'){
		qMap.plot1(bdm, data, qMap.levels = qMap.levels, qMap.cex = qMap.cex, qMap.bg = qMap.bg, qtitle = qtitle)
		return()
	}

	# get labels
	if (!is.null(labels)) labels <- as.numeric(as.factor(labels))

	# get var names
	if (is.null(colnames(data))) {
		colnames(data) <- paste('V', formatC(seq(ncol(data)), width = 2, flag = '0'), sep = '.')
	}

	# join labels & data
	if (!is.null(labels)) data <- cbind(labels, data)

	# check number of vars
	if (ncol(data) > 20) {
		data <- data[, 1:20]
		cat('+++ WARNING: plotting first ', ncol(data) - !is.null(labels), ' columns !, \n', sep='')
	}

	# get mapping
	l <- c(1, 2) + (layer -1) *2
	Y <- bdm$ptsne$Y[ , l]

	# set graphic environment (Att!! with this)
	layout.mtx <- t(cbind(layout.get(ncol(data)), layout.get(ncol(data))))
	layout.mtx[ , ] <- seq(length(layout.mtx))
	layout.mtx <- t(layout.mtx)
	layout(layout.mtx, widths = rep(2 /ncol(layout.mtx) *c(0.73, 0.27), ncol(layout.mtx)))
	c1 <- c(4.5, 3.5, 2.5, 1.5, 0.5)
	c3 <- c(8, 7, 6, 5, 4)
	parbdm.set(oma = c(c1[nrow(layout.mtx)], 1, c3[nrow(layout.mtx)], 1))

	# legend palette
	hmap.pltt <- c(pltt.heat(qMap.levels), '#DDDDDDFF')

	nulL <- lapply(seq(ncol(data)), function(j){
		if (j == 1 & !is.null(labels)) {
			X <- data[, j]
			if (is.null(class.pltt))
				pltt <- c(pltt.get(s = length(unique(X))), '#DDDDDDFF')
			else
				pltt <- class.pltt
		}
		else {
			# factor data
			X <- get.lvls(data[, j], qMap.levels)
			Q <- quantile(data[, j][!is.na(data[, j])], seq(0, 1, length.out = qMap.levels+1))
			pltt <- hmap.pltt
		}
		# plot q-maps
		par(mar = c(1.0, 1.0, 0.8, 0.4))
		if (!is.null(subset)) {
			# plot shadow
			plot(Y[-subset, ], xaxt = 'n', xlab = '', yaxt = 'n', ylab = '', xlim = range(Y[,1]), ylim = range(Y[,2]), col = pltt[length(pltt)], cex = qMap.cex, pch = 20, asp = 1, main = colnames(data)[j], cex.main = 0.6)
			# plot subset q-map
			points(Y[subset, ], col = pltt[X[subset]], cex = qMap.cex, pch = 20, asp = 1)
		}
		else {
			# plot q-map
			plot(Y, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = range(Y[, 1]), ylim = range(Y[, 2]), col = pltt[X], cex = qMap.cex, pch = 20, asp = 1, main = colnames(data)[j], cex.main = 0.6)
		}
		# plot title
		# title(colnames(data)[j], line = 1.0, outer = T, cex.sub = 0.8)
		# mtext(colnames(data)[j], side = 3, line = 2.0, at = 1.1 *min(Y[, 1]), cex = 1.0, pos = 'right')
		# text(0.8 *min(Y[, 1]), 0.9 *max(Y[, 2]), labels = colnames(data)[j], cex = 1.0)
		# plot legend
		par(mar = c(1, 0.1, 0.5, 0.1))
		plot(1, 1, xlab = '', ylab = '', xaxt = "n", yaxt = "n", bty = "n", type = "n")
		if (j == 1 & !is.null(labels)) {
			s <- length(unique(labels))
			lgnd.lbls <- formatC(seq(s), width = 3)
			legend('center', legend = lgnd.lbls[s:1], bty = 'n', pch = 15, cex = 0.6, pt.cex = 1.6, y.intersp = 0.7, col = pltt[s:1])
		} else {
			lgnd.lbls <- sapply(seq(qMap.levels), function(l) {
				if (l == 1 || Q[l] != Q[l-1]) {
					formatC(Q[l], format = 'e', digits = 2)
				} else {
					' '
				}
			})
			legend('center', legend = lgnd.lbls[qMap.levels:1], bty = 'n', pch = 15, cex = 0.6, pt.cex = 1.6, y.intersp = 0.7, col = pltt[qMap.levels:1])
		}
	})

	# fill layout
	if (length(ncol(data)) < max(layout.get(length(ncol(data))))) plot.null()
	# layout title
	if (is.null(qtitle)) qtitle <- bdm$dSet
	title(qtitle, outer = T, cex.main = 1.0)
	# reset graphic environment
	parbdm.def()
}

# ------------------------------------------------------------------------------
# +++ Precision map (quantile map of betas)
# ------------------------------------------------------------------------------

bdm_.pMap <- function(m, qMap.levels = 8, qMap.cex = 0.1, bg = '#000000')
{
	# get mapping
	Y <- m$ptsne$Y[, 1:2]

	# legend palette
	if (bg == '#000000')
		pltt <- c(pltt.heat(qMap.levels), '#DDDDDD')
	else
		pltt <- c(pltt.heat(qMap.levels), '#000000')

	#factor data
	data <- m$ppx$B[, 1]
	X <- get.lvls(data, qMap.levels)
	Q <- quantile(data, seq(0, 1, length.out = qMap.levels+1))

	# set graphic environment
	parbdm.set(oma = c(1.0, 1.0, 1.0, 1.0), mar = c(1.0, 1.0, 0.2, 0.4), bg = bg)
	layout(matrix(seq(2), nrow = 1), widths = c(0.73, 0.27))

	# plot q-map
	plot(Y, xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab = '', xlim = range(Y[, 1]), ylim = range(Y[, 2]), col = pltt[X], cex = qMap.cex, pch = 20, asp = 1)
	# plot legend
	par(mar = c(1.0, 0.1, 0.2, 0.1))
	plot(1, 1, xlab = '', ylab = '', xaxt = "n", yaxt = "n", bty = "n", type = "n")
	lgnd.lbls <- sapply(seq(qMap.levels), function(l) {
		if (l == 1 || Q[l] != Q[l-1]) {
			formatC(Q[l], format = 'e', digits = 2)
		} else {
			' '
		}
	})
	legend('center', legend = lgnd.lbls[qMap.levels:1], bty = 'n', pch = 15, cex = 1.2, pt.cex = 1.2, y.intersp = 1.1, col = pltt[qMap.levels:1], text.col = '#FFFFFF')

	# reset graphic environment
	parbdm.def()
}

qMap.plot1 <- function(bdm, data, qMap.levels = 8, qMap.cex = 0.3, qMap.bg = '#FFFFFF', qtitle =NULL)
{
	# get data
	if (is.null(data))
		return(message('+++ Error: no data given !'))
	# get mapping
	Y <- bdm$ptsne$Y[ , 1:2]
	# set graphic environment (Att!! with this)
	layout(matrix(1:2, nrow = 1), widths = c(0.73, 0.27))
	parbdm.set(oma = c(0.5, 0.5, 2.0, 0.5))
	# legend palette
	hmap.pltt <- c(pltt.heat(qMap.levels), '#DDDDDDFF')
	# factor data
	X <- get.lvls(data, qMap.levels)
	Q <- quantile(data[!is.na(data)], seq(0, 1, length.out = qMap.levels+1))
	pltt <- hmap.pltt
	par(mar = c(1.0, 1.0, 2.0, 0.4))
	# plot q-map
	plot(Y, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = range(Y[, 1]), ylim = range(Y[, 2]), col = pltt[X], cex = qMap.cex, pch = 20, asp = 1, main = qtitle, cex.main = 1.0)
	# plot legend
	par(mar = c(1, 0.1, 2.0, 0.1))
	plot(1, 1, xlab = '', ylab = '', xaxt = "n", yaxt = "n", bty = "n", type = "n")
	lgnd.lbls <- sapply(seq(qMap.levels), function(l) {
		if (l == 1 || Q[l] != Q[l-1]) {
			formatC(Q[l], format = 'e', digits = 2)
		} else {
			' '
		}
	})
	legend('center', legend = lgnd.lbls[qMap.levels:1], bty = 'n', pch = 15, cex = 1.0, pt.cex = 2.0, y.intersp = 1.0, col = pltt[qMap.levels:1])
	# reset graphic environment
	parbdm.def()
}
