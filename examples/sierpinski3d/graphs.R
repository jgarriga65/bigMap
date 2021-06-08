# -----------------------------------------------------------------------------
# +++ heat palette (this is a copy of the color.jet() function)
# -----------------------------------------------------------------------------

pltt.heat <- function(n, alpha = 1)
{
	if (length(n) > 1 || !is.finite(n))
		stop("'n' must be an integer positive value.", call. = FALSE)
	if (n < 1)
		stop("'n' must be an integer positive value.", call. = FALSE)
	if (length(alpha) > 1 || !is.finite(alpha))
		stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
	if (alpha < 0 || alpha > 1)
		stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
	alpha = round(255 * alpha)
	ramp = colorRamp(c("#00008F", "#00009F", "#0000AF", "#0000BF",
		"#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF",
		"#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF",
		"#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF",
		"#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF",
		"#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF",
		"#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60",
		"#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10",
		"#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00",
		"#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000",
		"#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000",
		"#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000",
		"#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000"),
		space = "Lab")
	rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
}

# ------------------------------------------------------------------------------
# +++ graph plot
# ------------------------------------------------------------------------------

graph.plot <- function(g, edges, cex = 0.8, lwd = 0.2, lvls = 8, bg = "#FFFFFF")
{
	Y <- g$ptsne$Y[, 1:2]
	par(oma = c(0.1, 0.1, 0.1, 0.1), mar = c(0, 0, 0, 0), bg = bg)
	plot(Y, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n', col='#888888', cex = cex, pch = 20, asp = 1.0)
	e.lens <- sqrt(apply(edges, 1, function(e) sum((Y[e[1], ] - Y[e[2], ])**2)))
	e.qntl <- quantile(e.lens, seq(0, 1, length.out = lvls+1))
	e.brks <- c(e.qntl[which(diff(e.qntl) != 0)], max(e.lens))
	e.lbls <- as.numeric(
		cut(e.lens, breaks=e.brks, labels=which(diff(e.qntl) != 0))
	)
	pltt <- pltt.heat(lvls)[lvls: 1]
	nulL <- sapply(sort(e.lens, decreasing = T, index.return = T)$ix, function(e)
	{
		i <- edges[e, 1]
		j <- edges[e, 2]
		lines(x = Y[c(i, j), 1], y = Y[c(i, j), 2], lwd=lwd, col=pltt[e.lbls[e]])
	})
}
