# ------------------------------------------------------------------------------
# +++ primes palette
# ------------------------------------------------------------------------------

primes.hpltt <- function(X, layer, offSet)
{
	# max. saturation value
	sat.depth <- max(X[, seq(2, ncol(X), by = 2)])
	# layer column
	col <- (layer -1) *2 +1
	# hue : prime factor
	hue.seq <- ((X[, col] +offSet) %%40) /45
	hue.seq[which(X[, col] == 0)] <- ((42.5 +offSet) %%40) /45
	# saturation: the highest the exponent the whiter the color
	sat.seq <- 1.0 -X[, (col +1)] /(2 *sat.depth)
	# transparency: zero exponent yields transparent color
	alpha.seq <- sapply(X[, (col +1)], function(w) ifelse(w == 0, '22', 'FF'))
	# palette
	paste(hsv(h = hue.seq, s = sat.seq, v = 0.9), alpha.seq, sep = '')
}

# ------------------------------------------------------------------------------
# +++ primes plot
# ------------------------------------------------------------------------------

primes.plot <- function(X, m, n.primes = 15, size.min = 2.0, size.max = 4.0, layers = 1:2, cex = 0.05, offSet = 28, noP = list(c(4, 8, 16, 32), c(9, 27, 81), c(25, 125, 625), c(49, 343)), bg = '#FFFFFF')
{
	Y <- m$ptsne$Y[, 1:2]
	# blurring factor to separate layers of different factors (though it is hardly visible)
	xBlurr <- sqrt(diff(range(Y[, 1]))**2 +diff(range(Y[, 2]))**2) /sqrt((724)**2 +(602)**2)
	xBlurr <- xBlurr /8
	# transparency for 3 layers
	alpha <- c('FF', '33', '11')

	par(mar = c(.2, .2, .2, .2), oma = c(.2, .2, .2, .2), bg = bg)

	plot(0, 0, xlab = '', ylab = '', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = range(Y[, 1]), ylim = range(Y[, 2]), col = '#FFFFFF')

	nulL <- sapply(layers, function(l) {
		# palette
		pltt <- primes.hpltt(X, l, offSet)
		# blurring factor
		rndN <- runif(nrow(X)) *ifelse((l == 1 || length(layers) == 1), 0, xBlurr) +1
		points(Y[, 1] *rndN, Y[, 2] *rndN, col = pltt, pch = 20, cex = cex)
	})

	if (n.primes > 0)
	{
		P <- primes::generate_n_primes(n.primes)
		if (n.primes > 15) {
			primes.seq <- P[seq(1, n.primes, length.out = 15)]
		} else {
			primes.seq <- P[seq(n.primes)]
		}

		primes.size <- seq(size.min, size.max, length.out = 16)[16: 1]
		text(x = m$ptsne$Y[1, 1], y = m$ptsne$Y[1, 2], labels = 0, cex = primes.size[1])
		nulL <- sapply(seq_along(primes.seq), function(k) {
			i <- primes.seq[k]
			# prime label
			text(x = m$ptsne$Y[i, 1], y = m$ptsne$Y[i, 2], labels = primes.seq[k], cex = primes.size[k])
			# prime lines
			if (k < 15) {
				j <- primes.seq[k +1]		# set for raw/P1G
				lines(x = m$ptsne$Y[c(i, j), 1], y = m$ptsne$Y[c(i, j), 2], lw = 0.8, col = 1)
			}
		})
	}

	if (!is.null(noP))
	{
		no.primes.size <- seq(size.min /3, size.min, length.out = 4)[4: 1]
		nulL <- sapply(seq_along(noP), function(k) {
			sapply(seq_along(noP[[k]]), function(ki) {
				i <- noP[[k]][[ki]]
				text(x = m$ptsne$Y[i, 1], y = m$ptsne$Y[i, 2], labels = i, cex = no.primes.size[ki], col = 2)
			})
		})
	}
}
