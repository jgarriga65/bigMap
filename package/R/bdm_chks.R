
# -----------------------------------------------------------------------------
# +++ nSS check functions
# -----------------------------------------------------------------------------

bdm.comp.betas <- function(B1, B2)
{
	layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = T))
	y.lim <- c(min(c(B1, B2)), max(c(B1, B2)))
	plot(seq(length(B1)), sort(B1), type = 'l', col = 4, lwd = 2, ylim = y.lim, xlab = 'i', ylab = 'Beta_i')
	lines(seq(length(B1)), sort(B2), col = 2, lwd = 2)
	legend('topleft', legend = c(substitute(B1), substitute(B2)), col = c(4, 2), bty = 'n', lty = 1)
	B4 <- (B1 -max(B1))/diff(range(B1))
	B5 <- (B2 -max(B2))/diff(range(B2))
	y.lim <- c(min(c(B4, B5)), max(c(B4, B5)))
	plot(seq(length(B1)), sort(B4), type = 'l', col = 4, lwd = 2, ylim = y.lim, xlab = 'i', ylab = 'Beta_i')
	lines(seq(length(B1)), sort(B5), col = 2, lwd = 2)
	B7 <- (B1 -mean(B1)) - (B2 -mean(B2))
	plot(seq(length(B1)), B7, type = 'l', col = 4, lwd = 1.0, xlab = 'i', ylab = 'B1-B2')
}

bdm.chk.nnSize <- function(D, m.ppx, is.distance, is.sparse, normalize, threads, layers)
{
	zSize <- round(nrow(D) *layers /threads, 0)
	nnSize <- ptsne.nnSize(m.ppx$ppx, m.ppx$xppx, layers, threads, zSize)
	cat('+++ zSize ', zSize, ', expected nnSize ', nnSize, ', nnRatio', round(nnSize/zSize, 4), '\n')
	Xbm <- as.big.matrix(t(D), type = 'double')
	if (normalize) centerScale(Xbm@address, is.distance, is.sparse)
	Bbm <- as.big.matrix(t(m.ppx$B), type = 'double')
	thread.nnSize <- as.numeric(sapply(seq(10), function(round) {
		zIdx <- sample(seq(ncol(Xbm)))[1: zSize]
		nnSS <- nnSS_chk(Xbm@address, Bbm@address, zIdx, is.distance, is.sparse, nnSize)
		cat('... max.', max(nnSS), ', at', zIdx[which.max(nnSS)], ', Beta.', Bbm[1, zIdx[which.max(nnSS)]], 'Lower.', length(which(nnSS < nnSize)), 'Zero.', length(which(nnSS < 3)), '\n')
		nnSS
	}))
	print(summary(thread.nnSize))
	thread.nnSize
}

thread.chk.nnSize <- function()
{
	if (thread.rank != 0) {
		# Att!! X transposed
		zIdx <- sample(seq(ncol(Xbm)))[1:zSize]
		nnSS_chk(Xbm@address, Bbm@address, zIdx, is.distance, is.sparse, nnSize)
	}
}

bdm.chk.aff <- function(D, chkB, is.distance, is.sparse, zSize, threads = 100, layers = 2)
{
	nnSize <- bdm.nnSize(chkB$ppx, chkB$xppx, layers, threads, zSize)
	cat('+++ zSize ', zSize, ', expected nnSize ', nnSize, ', nnRatio', round(nnSize/zSize, 4), '\n')
	X <- as.big.matrix(t(D), type = 'double')
	B <- as.big.matrix(t(chkB$B), type = 'double')
	zIdx <- sample(seq(3:nrow(D)))[1: zSize]
	zIdx[1] <- 1
	zAff <- aff_chk(X@address, B@address, zIdx, is.distance, is.sparse, nnSize)
	list(rows = zIdx, mtx = zAff)
}

bdm.chk.time <- function(D, chkB, is.distance, is.sparse, normalize, threads, layers, useEx = 1.0, theta = .0, epochs = 3, iters = NULL)
{
	m <- list()
	m$dSet <- 'chkTimes'
	m$is.distance <- is.distance
	m$is.sparse <- is.sparse
	m$normalize <- normalize
	#
	m$Xdata$whiten <- 0
	m$Xdata$input.dim <- NULL
	#
	m$Xbeta$ppx <- chkB$ppx
	m$Xbeta$nnSize <- chkB$xppx
	m$Xbeta$B <- chkB$B
	#
	m$ptsne$threads <- threads
	m$ptsne$layers <- layers
	m$ptsne$rounds <- 1
	m$ptsne$theta <- theta
	m$ptsne$alpha <- 0.5
	m$ptsne$useEx <- useEx
	#
	bdm.restart(D, m, epochs = epochs, iters = iters)
}

# checks that Y.in is chunked correctly and Y.out is pooled correctly
bdm.chkmpi.simple <- function(nX, threads, layers)
{
	chnk.brks <- round(seq(1, nX +1, length.out = (threads +1)), 0)
	thrd.size <- round(nX /threads *layers, 0)
	thrd.brks <- lapply(seq(threads), function(z)
	{
		t(sapply(seq(layers), function(l) {
			a <- z + (l-1)
			if (a > threads) a <- a %% threads
			# Att!! C++ indexes
			c(chnk.brks[a], chnk.brks[a+1]) -1
		}))
	})
	Z.list <- lapply(seq(threads), function(z) {
		matrix(0, sum(apply(thrd.brks[[z]], 1, diff)), 3)
	})
	# initial mapping (by default, random embedding on a unity disk)
	Y.ini <- ptsne.init(nX, layers)
	Y.out <- ptsne.init(nX, layers) *10
	# resample dataset. Att!! C++ indexes
	I <- sample(seq(nX)) -1
	cl <- cluster.start(threads -1, 'SOCK')
	clusterEvalQ(cl, library(Rcpp))
	clusterEvalQ(cl, sourceCpp('~/bigMap/bigMap_4.5.0/src/zTSNE.cpp'))
	# special initialization for thread.rank != 0
	clusterEvalQ(cl,
		if (thread.rank != -1) {
			w <- new.env()
			w$zI <- numeric()
			w$zY <- numeric()
		})
	zChnks(Z.list, Y.ini, I, thrd.brks)
	nulL <- clusterApply(cl, Z.list, chk.ztsne)
	# pool partial mappings from workers
	zMap.list <- clusterEvalQ(cl, if (thread.rank != -1) w$zY)
	# restructure global mapping
	updateY(Y.out, I, zMap.list, thrd.brks)
	cluster.stop(cl)
	return(list(In = Y.ini, Out = Y.out))
}

chk.ztsne <- function(zChnk)
{
	if (thread.rank != -1) {
		w$zI <- zChnk[, 1]
		w$zY <- zChnk[, 2:3]
		chk_zTSNE(w$zY, w$zI)
	}
}

bdm.chk.rowDist <- function(X, is.distance, is.sparse, row, ppx)
{
	Xbm <- as.big.matrix(t(X), type = 'double')
	centerScale(Xbm@address)
	rowDist(Xbm@address, is.distance, is.sparse, (row -1), ppx)
}

bdm.chk.rowBeta <- function(X, is.distance, is.sparse, row, ppx, xppx)
{
	Xbm <- as.big.matrix(t(X), type = 'double')
	centerScale(Xbm@address)
	rowBeta(Xbm@address, is.distance, is.sparse, (row -1), ppx, xppx)
}

bdm.chk.efc <- function(X, m, threads, thread.rank = 1, iters = 10)
{
	bigX <- as.big.matrix(t(X), type = 'double')
	centerScale(bigX@address, m$is.distance, m$is.sparse)

	bigB <- as.big.matrix(t(m$Xbeta$B), type = 'double')
	bigY <- as.big.matrix(t(m$ptsne$Y[, 1:2]), type = 'double')

	eRange <- c(range(m$ptsne$Y[, 1]), range(m$ptsne$Y[, 2]))
	nnSize <- m$Xbeta$ppx *m$Xbeta$nnSize
	EFC(thread.rank, threads, bigX@address, m$is.distance, m$is.sparse, bigB@address, bigY@address, iters)
}
