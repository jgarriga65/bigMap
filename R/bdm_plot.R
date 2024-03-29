# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# +++ ptSNE plots
# -----------------------------------------------------------------------------

#' ptSNE cost & size plot.
#'
#' @param bdm A \var{bdm} data mapping instance, or a list of them to make a comparative plot.
#'
#' @param x.lim X-axis (epochs) limits. Default value is \var{x.lim = NULL}.
#'
#' @return None.
#'
#' @examples
#'
#' bdm.example()
#' bdm.cost(ex$map)

bdm.cost <- function(bdm, x.lim = NULL)
{
	if (is.null(names(bdm))) bdm.list <- bdm
	else bdm.list <- list(bdm)

	# set graphic environment
	parbdm.set(mar=c(4.5,4.5,3,4.5), mgp=c(1.0,0.6,0), cex.axis=1.2)
	layout(layout.get(length(bdm.list)))

	nulL <- lapply(bdm.list, function(bdm)
	{
		ptsne.cost(bdm, x.lim = x.lim)
	})

	# fill layout
	if (length(bdm.list) < max(layout.get(length(bdm.list)))) plot.null()
	# reset graphic environment
	parbdm.def()

}

# +++ Plot cost/size internal function.

ptsne.cost <- function(bdm, x.lim = NULL, movie = F, mtext.cex = 1.2)
{
	if (!is.null(bdm$ptsne$Y))
	{
		if (bdm$ptsne$threads == 1) bdm$ptsne$cost <- t(bdm$ptsne$cost)
		cost.clr <- brewer.pal(5, 'Purples')
		if (is.null(x.lim)) x <- 0:(ncol(bdm$ptsne$cost)-1)
		else x <- (x.lim[1]: x.lim[2]) -1
		y <- x +1
		# cost by threads
		if (bdm$ptsne$threads == 1) {
			plot(x, bdm$ptsne$cost[1, y], type='l', col=cost.clr[5], axes=F, xlab='', ylab='', ylim = range(bdm$ptsne$cost[1, y]), lwd = 1.5)
		} else {
			thread <- 1
			plot(x, bdm$ptsne$cost[thread, y], type='l', col=cost.clr[2], axes=F, xlab='', ylab='', ylim=range(bdm$ptsne$cost))
			nulL <- sapply(seq(2, nrow(bdm$ptsne$cost)), function(thread)
			{
				lines(x, bdm$ptsne$cost[thread, y], col=cost.clr[2])
			})
			# average cost
			lines(x, apply(bdm$ptsne$cost[ ,y], 2, mean), col=cost.clr[5], lwd=1.5)
		}
		axis(side=1, at=pretty(range(x)), tick=T)
		if (!movie) mtext("epochs", side=1, line=ifelse((mtext.cex==1.2), 3, 2), cex=mtext.cex)
		axis(side=2, at=pretty(range(bdm$ptsne$cost[1, y])), tick=T, las=1, col=cost.clr[5])
		mtext("Cost", side=2, line=ifelse((mtext.cex==1.2), 3, 2), col=cost.clr[5], cex=mtext.cex)
		# layers' average size
		if (!is.null(bdm$ptsne$size))
		{
			par(new = T)
			size.clr <- brewer.pal(5, 'PuRd')
			plot(x, bdm$ptsne$size[y], type = 'l', col = size.clr[5], lwd = 1.5, axes = F, xlab = '', ylab= '', ylim = range(bdm$ptsne$size[y]))
			axis(side=4, at=pretty(range(bdm$ptsne$size[y])), tick=T, las=1, col=size.clr[5])
			mtext("Size", side=4, line=ifelse((mtext.cex==1.2), 3, 1.5), col=size.clr[5], cex=mtext.cex)
		}
	}
}

#' Plot ptSNE (low-dimensional embedding)
#'
#' @param bdm A \var{bdm} data mapping instance, or a list of them to make a comparative plot.
#'
#' @param ptsne.cex The size of the mapped data-points in the ptSNE plot. Default value is \var{ptsne.cex = 0.5}.
#'
#' @param ptsne.bg The background colour of the ptSNE plot. Default value is \var{ptsne.bg = #FFFFFF} (white).
#'
#' @param class.lbls A vector of class labels to show in the ptSNE plot. Default value is \var{class.lbls = NULL}. If \var{is.null(class.lbls)} and \var{!is.null(bdm$lbls)} the later will be used. If \var{!is.null(bdm$lbls)} and \var{!is.null(bdm$wtt)} (i.e. there is a clustering) cluster labels are used by default.
#'
#' @param class.pltt A colour palette to show class labels in the ptSNE plot. If \var{!is.null(bdm$wtt)} cluster labels are used by default, else if \var{!is.null(bdm$lbls)} are used by default. If \var{ptsne.pltt = NULL} (default value) the default palette is used.
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @return None.
#'
#' @examples
#'
#' bdm.example()
#' bdm.ptsne.plot(ex$map)

bdm.ptsne.plot <- function(bdm, ptsne.cex = 0.5, ptsne.bg = '#FFFFFF', class.lbls = NULL, class.pltt = NULL, layer = 1)
{
	if (is.null(names(bdm))) bdm.list <- bdm
	else bdm.list <- list(bdm)
	# set graphic environment
	parbdm.set(oma = c(0.8, 0.8, 0.8, 0.8), mar = c(2.8, 2.8, 0.5, 0.5), mgp=c(1.8,0.6,0), cex.axis=1.0)
	layout(layout.get(length(bdm.list)))
	nulL <- lapply(bdm.list, function(bdm)
	{
		if (!is.null(bdm$ptsne)) {
			ptsne.plot(bdm, layer = layer, cex = ptsne.cex, bg = ptsne.bg, lbls = class.lbls, pltt = class.pltt)
		}
		else {
			plot.null()
			return(message('+++ Error: no ptSNE found \n'))
		}
	})
	# fill layout
	if (length(bdm.list) < max(layout.get(length(bdm.list)))) plot.null()
	# reset graphic environment
	parbdm.def()
}

# ------------------------------------------------------------------------------
# +++ ptSNE scatterplot (internal)
# ------------------------------------------------------------------------------

ptsne.plot <- function(bdm, lbls = NULL, pltt = NULL, cex = 0.3, bg = '#FFFFFF', layer = 1)
{
	if (!is.null(lbls)) {
		L <- lbls
	}
	else if (!is.null(bdm$lbls)) {
		L <- bdm$lbls
	}
	else {
		L <- bdm.labels(bdm, layer = layer)
	}
	if (is.null(pltt)) {
		if (!is.null(bdm$merge)) {
			pltt <- rep(0, bdm$wtt[[layer]]$s)
			pltt[unique(L)] <- pltt.get(length(unique(L)))
		} else {
			pltt <- pltt.get(length(unique(L)))
		}
	}
	if (layer == -1) {
		Y1 <- apply(bdm$ptsne$Y[, seq(1, bdm$ptsne$layers *2, by = 2)], 1, mean)
		Y2 <- apply(bdm$ptsne$Y[, seq(2, bdm$ptsne$layers *2, by = 2)], 1, mean)
	}
	else {
		Y1 <- bdm$ptsne$Y[, (1 +(layer -1) *2)]
		Y2 <- bdm$ptsne$Y[, (2 +(layer -1) *2)]
	}

	par(bg = bg)
	plot(Y1, Y2, xlab = 'Y1', ylab = 'Y2', col = pltt[L], pch = 20, cex = cex, cex.lab = 1.0)

}

#' Plot paKDE (density landscape)
#'
#' @param bdm A \var{bdm} data mapping instance, or a list of them to make a comparative plot.
#'
#' @param pakde.pltt A colour palette to show levels in the paKDE plot. By default (\var{pakde.pltt = NULL}) the default palette is used.
#'
#' @param pakde.lvls The number of levels of the density heat-map (16 by default).
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @return None.
#'
#' @examples
#'
#' bdm.example()
#' bdm.pakde.plot(ex$map)

bdm.pakde.plot <- function(bdm, pakde.pltt = NULL, pakde.lvls = 16, layer = 1)
{
	if (is.null(names(bdm))) bdm.list <- bdm
	else bdm.list <- list(bdm)
	# set graphic environment
	parbdm.set(oma = c(0.8, 0.8, 0.8, 0.8), mar = c(2.8, 2.8, 0.5, 0.5), mgp=c(1.8,0.6,0), cex.axis=1.0)
	layout(layout.get(length(bdm.list)))
	nulL <- lapply(bdm.list, function(bdm)
	{
		if (is.null(bdm$pakde)) {
			plot.null()
			return(message('+++ Error: no ptSNE found \n'))
		}
		else {
			pakde <- bdm$pakde[[layer]]
			plot.pakde(pakde, pakde.pltt, pakde.lvls)
		}
	})
	# fill layout
	if (length(bdm.list) < max(layout.get(length(bdm.list)))) plot.null()
	# reset graphic environment
	parbdm.def()
}

#' Plot WTT (clustering)
#'
#' @param bdm A \var{bdm} data mapping instance, or a list of them to make a comparative plot.
#'
#' @param pakde.pltt A colour palette to show levels in the paKDE plot. By default (\var{pakde.pltt = NULL}) the default palette is used.
#'
#' @param pakde.lvls The number of levels of the density heat-map (16 by default).
#'
#' @param wtt.lwd The width of the watertrack lines (as set in \var{par()}).
#'
#' @param plot.peaks Default value is \var{plot.peaks = TRUE}. If TRUE the density peaks are indicated.
#'
#' @param labels.cex If \var{plot.peaks} is TRUE, the size of the peaks labels (as set in \var{par()}). Default value is \var{labels.cex = 1.0}.
#'
#' @param layer The ptSNE output layer. Default value is \var{layer = 1}.
#'
#' @return None.
#'
#' @examples
#'
#' bdm.example()
#' bdm.wtt.plot(ex$map)

bdm.wtt.plot <- function(bdm, pakde.pltt = NULL, pakde.lvls = 16, wtt.lwd = 1.0, plot.peaks = T, labels.cex = 1.0, layer = 1)
{
	if (is.null(names(bdm))) bdm.list <- bdm
	else bdm.list <- list(bdm)
	# set graphic environment
	parbdm.set(oma = c(0.8, 0.8, 0.8, 0.8), mar = c(2.8, 2.8, 0.5, 0.5), mgp=c(1.8,0.6,0), cex.axis=1.0)
	layout(layout.get(length(bdm.list)))
	nulL <- lapply(bdm.list, function(bdm)
	{
		if (is.null(bdm$pakde)) {
			plot.null()
			return(message('+++ Error: no ptSNE found \n'))
		}
		else {
			pakde <- bdm$pakde[[layer]]
			plot.pakde(pakde, pakde.pltt, pakde.lvls)
			if (!is.null(bdm$wtt) && layer <= length(bdm$wtt))
			{
				wtt <- bdm$wtt[[layer]]
				if (!is.null(bdm$merge)) {
					plot.wtt(pakde, bdm$merge$C, wtt$grid, 2*wtt.lwd, '#555555')
					wtt.lwd <- wtt.lwd * 0.5
				}
				plot.wtt(pakde, wtt$C, wtt$grid, wtt.lwd, '#CCCCCC')
				if (plot.peaks) {
					if (!is.null(bdm$merge)) C <- bdm$merge$C
					else C <- wtt$C
					wtt.peaks(pakde, wtt, C, labels.cex)
				}
			}
		}
	})
	# fill layout
	if (length(bdm.list) < max(layout.get(length(bdm.list)))) plot.null()
	# reset graphic environment
	parbdm.def()
}

# ------------------------------------------------------------------------------
# +++ plot pakde (internal)
# ------------------------------------------------------------------------------

plot.pakde <- function(pakde, pltt, lvls)
{
	if (is.null(pltt)) pltt <- pltt.pakde(lvls)
	image(pakde$x, pakde$y, pakde$z, col = pltt, xaxt='n', yaxt='n', xlab='', ylab='', )
}

# ------------------------------------------------------------------------------
# +++ plot wtt.lines (internal)
# ------------------------------------------------------------------------------

plot.wtt <- function(pakde, C, grid, lwd, col)
{
	nulL <- sapply(seq_along(C), function(n)
	{
		n2c <- as.numeric(grid_n2cell(n-1, grid)) +1
		n.cross <- as.numeric(grid_cross(n-1, grid)) + 1
		nulL <- sapply(n.cross, function(m)
		{
			if (m > n && C[n] != C[m]) {
				m2c <- as.numeric(grid_n2cell(m-1, grid)) +1
				if (n2c[1] != m2c[1]) {
					lines(pakde$x[c(m2c[1], m2c[1])], pakde$y[c(n2c[2], (n2c[2]+1))], col=col, lwd=lwd)
				}
				if (n2c[2] != m2c[2]) {
					lines(pakde$x[c(n2c[1], (n2c[1]+1))], pakde$y[c(m2c[2], m2c[2])], col=col, lwd=lwd)
				}
			}
		})
	})
}


# ------------------------------------------------------------------------------
# +++ add peaks to wtt.lines plot (internal)
# ------------------------------------------------------------------------------

wtt.peaks <- function(pakde, wtt, C, labels.cex)
{
	peaks <- unique(C)
	points(pakde$x[wtt$M[peaks, 1]], pakde$y[wtt$M[peaks, 2]],  col='#000000FF', cex=1, pch=17)
	if (labels.cex > 0) {
		text(pakde$x[wtt$M[peaks, 1]], pakde$y[wtt$M[peaks, 2]], labels = peaks, pos = 3, cex = labels.cex)
	}
}
