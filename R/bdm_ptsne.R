# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ptsne.scheme <- function(nX, layers, threads, ppx, xppx = 3.0, base_ppx = 0.03)
{
	zSize <- nX *layers /threads
	base_nnS <- ptsne.nnSize(base_ppx *nX, xppx, layers, threads, zSize)
	nnS <- round(ptsne.nnSize(ppx, xppx, layers, threads, zSize), 0)
	# log2(2**20) = log2(1048576) = 20; 2**7 = 128
	# epochs <- floor(128 *log2(nX) /20 *log2(base_nnS) /log2(nnS))
	# epochs <- floor(64 *log2(nX) /20)
	epochs <- floor(4 *log(nX))
	# log2(2**14) = log2(16384) = 14
	# iters <- ceiling(32 *log2(zSize) /14)
	iters <- ceiling(4 *log(zSize))
	scheme.list <- list(epochs = epochs, iters = iters)
	# if (ppx >= (base_ppx *nX)) {
	# 	scheme.list <- list(epochs = epochs, iters = iters)
	# }
	# else {
	# 	scheme.list <- list(epochs = iters, iters = epochs)
	# }
	return(scheme.list)
}

ptsne.nnSize <- function(ppx, xppx, layers, threads, zSize)
{
	if (threads /layers == 1) nnSize <- min(3 *ppx, zSize -1)
	else nnSize <- round(min(xppx *ppx *layers /threads, zSize -1), 0)
	return(nnSize)
}

thread.chk.nnSize <- function()
{
	if (thread.rank != 0) {
		# Att!! X transposed (C++ indexes)
		zIdx <- sample(seq(ncol(Xbm)))[1:zSize] -1
		# chnk.brks <- round(seq(1, ncol(Xbm) +1, length.out = (threads +1)), 0)
		# if (thread.rank < threads)
		# 	zIdx <- chnk.brks[thread.rank]:(chnk.brks[thread.rank+2] -1)
		# else
		# 	zIdx <- c(chnk.brks[thread.rank]:(chnk.brks[thread.rank+1] -1), chnk.brks[1]:(chnk.brks[2] -1))
		nnSS_chk(Xbm@address, Bbm@address, zIdx, is.distance, is.sparse, nnSize)
	}
}

ptsne.get <- function(cl, bdm, info)
{
	bdm.list <- lapply(bdm$ppx, function(ppx) {
		m <- bdm
		m$ppx <- ppx
		# +++ export betas
		cat('+++ exporting betas, ppx ', m$ppx$ppx, '\n')
		m$t$exportBetas <- system.time(
			Xbeta.exp(cl, m$ppx$B)
		)
		print(summary(m$ppx$B[, 1]))
		print(m$t$exportBetas)
		# +++ check nnSize
		cat('+++ checking nnSize \n')
		m$t$nnSize <- system.time({
			threads <- m$ptsne$threads
			layers  <- m$ptsne$layers
			zSize <- round(m$nX /threads *layers, 0)
			nnSize <- ptsne.nnSize(m$ppx$ppx, m$ppx$xppx, layers, threads, zSize)
			clusterExport(cl, c('nnSize', 'layers', 'zSize'), envir=environment())
			chk.list <- clusterCall(cl, thread.chk.nnSize)
			sumP <- mean(unlist(lapply(chk.list, function(thrd) thrd$P)))
			cat('... sumP', sumP, '\n')
			m$ppx$nnSize <- summary(unlist(lapply(chk.list, function(thrd) thrd$nnS)))
			print(m$ppx$nnSize)
		})
		print(m$t$nnSize)
		# +++ perform t-SNE
		cat('+++ computing ptSNE \n')
		m$t$ptsne <- system.time({
			if (attr(cl[[1]], 'class') == 'SOCKnode') {
				m <- sckt.ptsne(cl, m, info)
			} else {
				m <- mpi.ptsne(cl, m, info)
			}
		})
		print(m$t$ptsne)
		m
	})
	return(bdm.list)
}

ptsne.restart <- function(cl, bdm, info)
{
	# +++ export betas
	cat('+++ exporting betas, ppx ', bdm$ppx$ppx, '\n')
	Xbeta.exp(cl, bdm$ppx$B)
	# +++ perform t-SNE
	cat('+++ computing ptSNE \n')
	t <- system.time({
		if (attr(cl[[1]], 'class') == 'SOCKnode') {
			bdm <- sckt.ptsne(cl, bdm, info)
		} else {
			bdm <- mpi.ptsne(cl, bdm, info)
		}
	})
	print(t)
	return(bdm)
}


# -----------------------------------------------------------------------------
# +++ distributed ptSNE implementation with shared memory
# -----------------------------------------------------------------------------

sckt.ptsne <- function(cl, bdm, progress)
{
	# setup parameters
	threads  <- bdm$ptsne$threads
	layers   <- bdm$ptsne$layers
	ppx      <- bdm$ppx$ppx
	theta    <- bdm$ptsne$theta
	momentum <- bdm$ptsne$momentum
	qDecay   <- bdm$ptsne$qDecay
	gain     <- bdm$ptsne$gain
	#
	chnk.brks <- round(seq(1, bdm$nX +1, length.out = (threads +1)), 0)
	zSize <- round(bdm$nX /threads *layers, 0)

	nnSize <- ptsne.nnSize(ppx, bdm$ppx$xppx, layers, threads, zSize)
	bdm$ptsne$nnSize <- nnSize

	scheme <- ptsne.scheme(bdm$nX, layers, threads, ppx, bdm$ppx$xppx)
	epochs <- scheme$epochs
	iters <- scheme$iters
	# the next condition is meant only for development purposes
	if (!is.null(bdm$epochs)) epochs <- bdm$epochs
	if (!is.null(bdm$iters)) iters <- bdm$iters
	bdm$ptsne$epochs <- epochs
	bdm$ptsne$iters <- iters

	# export setup parameters
	clusterExport(cl, c('layers', 'nnSize', 'zSize', 'theta', 'gain'), envir = environment())
	clusterExport(cl, c('progress'), envir = environment())

	# initial mapping (by default, random circular embedding of radius 1)
	if (!is.null(bdm$ptsne$Y)) {
		Ybm <- as.big.matrix(bdm$ptsne$Y, type='double')
	}
	else {
		Ybm <- as.big.matrix(ptsne.init(bdm$nX, layers), type='double')
	}
	Ybm.dsc <- describe(Ybm)

	# embedding size
	eSize <- rep(0, (epochs +1))
	eSize[1] <- sqrt(sum(apply(apply(Ybm[, 1:2], 2, range), 2, diff)**2))

	# learning Rate
	lRate <- 2.0 *(eSize[1] +1.0 /eSize[1]) *log(zSize *nnSize)
	# momentum
	alpha <- momentum
	clusterExport(cl, c('lRate', 'alpha'), envir = environment())

	# row sampling indexes (substract 1 to convert to C++ indexes !!)
	Ibm <- big.matrix(bdm$nX, 1, type='integer')
	Ibm.dsc <- describe(Ibm)
	# epoch/thread cost matrix
	Cbm <- as.big.matrix(matrix(0, threads, (epochs +1)), type='double')
	Cbm.dsc <- describe(Cbm)

	# attach bigmatrices to workers
	clusterExport(cl, c('Ybm.dsc'), envir=environment())
	clusterEvalQ(cl, Ybm <- attach.big.matrix(Ybm.dsc))
	clusterExport(cl, c('Ibm.dsc'), envir=environment())
	clusterEvalQ(cl, Ibm <- attach.big.matrix(Ibm.dsc))
	clusterExport(cl, c('Cbm.dsc'), envir=environment())
	clusterEvalQ(cl, Cbm <- attach.big.matrix(Cbm.dsc))

	# special initialization for thread.rank == 0
	clusterEvalQ(cl,
		if (thread.rank == 0) {
			w <- new.env()
			w$mapp.list <- list()
		})

	clusterExport(cl, c('epochs', 'iters'), envir = environment())

	# compute initial cost
	clusterEvalQ(cl, epoch <- 0)
	clusterEvalQ(cl, iters <- 0)
	# resample dataset (C++ indexes)
	Ibm[ ] <- sample(seq(0, (bdm$nX -1)))

	# perform ptSNE
	nulL <- clusterCall(cl, sckt.ztsne)

	# starting embedding size
	eSize[1] <- sqrt(sum(apply(apply(Ybm[, 1:2], 2, range), 2, diff)**2))

	# learning Rate
	lRate <- 2.0 *(eSize[1] +1.0 /eSize[1]) *log(zSize *nnSize)
	clusterExport(cl, c('lRate'), envir = environment())

	# reset number of iterations
	clusterExport(cl, c('epochs', 'iters'), envir=environment())

	# start
	t0 <- Sys.time()

	# report starting information
	avgCost <- mean(Cbm[, 1])
	if (progress >= 0) {
		nulL <- ptsne.info(threads, zSize, nnSize, epochs, iters, 0, avgCost, eSize[1], t0, theta)
	}

	for (e in seq(epochs)) {
		# epoch start-time
		te <- Sys.time()
		# resample dataset (C++ indexes)
		Ibm[ ] <- sample(seq(0, (bdm$nX -1)))
		# perform ptSNE
		nulL <- clusterCall(cl, sckt.ztsne)
		# security break control
		if (mean(Cbm[, (e +1)]) < 0) break
		# embedding size
		eSize[(e +1)] <- sqrt(sum(apply(apply(Ybm[, 1:2], 2, range), 2, diff)**2))
		# learning Rate
		lRate <- 2.0 *(eSize[(e +1)] +1.0 /eSize[(e +1)]) *log(zSize *nnSize)
		# momentum
		if (qDecay)
			alpha <- momentum *(1 -e /scheme$epochs)**2
		else
			alpha <- momentum *(1 -e /scheme$epochs)
		clusterExport(cl, c('lRate', 'alpha'), envir = environment())
		# report status
		if (progress >=0) {
			avgCost <- mean(Cbm[, (e +1)])
			epoch.info(e, epochs, avgCost, eSize[(e +1)], t0, te, lRate)
		}
	}

	bdm$ptsne$Y <- as.matrix(Ybm[ , ])
	bdm$ptsne$cost <- as.matrix(Cbm[, 1:(epochs +1)])
	bdm$ptsne$size <- eSize

	if (progress == 2) {
		bdm$progress <- clusterEvalQ(cl, if (thread.rank == 0) w$mapp.list)[[1]]
	}

	# report status
	avgCost <- mean(Cbm[, e +1])
	bdm$t[['epoch']] <- ptsne.info(threads, zSize, nnSize, epochs, iters, e, avgCost, eSize[e +1], t0, theta)

	return(bdm)

}

# -----------------------------------------------------------------------------
# +++ ptSNE SCKT thread functions
# -----------------------------------------------------------------------------

sckt.ztsne <- function()
{
	epoch <<- epoch +1	# Att.!! global variable
	if (thread.rank == 0) {
		# save current embedding to make movie
		if (progress == 2) {
			w$mapp.list[[epoch]] <- list(epoch = epoch, Y = as.matrix(Ybm[ , 1:2]))
		}
	}
	else {
		zCost <- sckt_zTSNE(thread.rank, epoch, threads, layers, Xbm@address, Bbm@address, Ybm@address, Ibm@address, iters, nnSize, theta, lRate, alpha, gain, is.distance, is.sparse)
		Cbm[thread.rank, epoch] <- zCost
	}
}

# -----------------------------------------------------------------------------
# +++ ptSNE MPI
# -----------------------------------------------------------------------------

mpi.ptsne <- function(cl, bdm, progress)
{
	# setup parameters
	ppx      <- bdm$ppx$ppx
	threads  <- bdm$ptsne$threads
	layers   <- bdm$ptsne$layers
	theta    <- bdm$ptsne$theta
	momentum <- bdm$ptsne$momentum
	qDecay   <- bdm$ptsne$qDecay
	gain     <- bdm$ptsne$gain
	#
	chnk.brks <- round(seq(1, bdm$nX +1, length.out = (threads +1)), 0)
	zSize <- round(bdm$nX /threads *layers, 0)

	nnSize <- ptsne.nnSize(ppx, bdm$ppx$xppx, layers, threads, zSize)
	bdm$ptsne$nnSize <- nnSize

	scheme <- ptsne.scheme(bdm$nX, layers, threads, ppx, bdm$ppx$xppx)
	epochs <- scheme$epochs
	iters <- scheme$iters
	# the next condition is meant only for development purposes
	if (!is.null(bdm$epochs)) epochs <- bdm$epochs
	if (!is.null(bdm$iters)) iters <- bdm$iters
	bdm$ptsne$epochs <- epochs
	bdm$ptsne$iters <- iters

	# export ptSNE setup parameters
	clusterExport(cl, c('layers', 'nnSize', 'zSize', 'theta', 'gain'), envir = environment())
	clusterExport(cl, c('progress'), envir = environment())

	# thread segments: row/col indexes in Y
	zBrks <- lapply(seq(threads), function(z)
	{
		t(sapply(seq(layers), function(l) {
			a <- z + (l-1)
			if (a > threads) a <- a %% threads
			# Att!! C++ indexes
			c(chnk.brks[a], chnk.brks[a+1]) -1
		}))
	})

	# initialize Z.list to map I and Y chunks to workers
	Z.list <- lapply(seq(threads), function(z) {
		matrix(0, sum(apply(zBrks[[z]], 1, diff)), 3)
	})

	# initial mapping (by default, random embedding on a unity disk)
	if (!is.null(bdm$ptsne$Y)) {
		Y <- bdm$ptsne$Y
	}
	else {
		Y <- ptsne.init(bdm$nX, layers)
	}

	# initialize cost&size functions
	eCost <- matrix(0, nrow = threads, ncol = epochs +1)
	eSize <- rep(0, (epochs +1))
	eSize[1] <- sqrt(sum(apply(apply(Y[, 1:2], 2, range), 2, diff)**2))

	# learning Rate
	lRate <- 2.0 *(eSize[1] +1.0 /eSize[1]) *log(zSize *nnSize)
	# momentum
	alpha <- momentum
	clusterExport(cl, c('lRate', 'alpha'), envir = environment())

	# special initialization for thread.rank != 0
	clusterEvalQ(cl,
		if (thread.rank != 0) {
			w <- new.env()
			w$zI <- numeric()
			w$zY <- numeric()
		})

	clusterExport(cl, c('epochs', 'iters'), envir = environment())

	# +++ compute initial cost&size
	clusterEvalQ(cl, epoch <- 0)
	clusterEvalQ(cl, iters <- 0)
	# resample dataset. Att!! C++ indexes
	I <- sample(seq(bdm$nX)) -1
	# get I&Y chunks
	zChnks(Z.list, Y, I, zBrks)

	# pool cost from workers
	# NULL returned by thread.rank == 0 is removed by unlist
	eCost[, 1] <- unlist(clusterApply(cl, Z.list, mpi.ztsne))

	# update lRate to workers
	eSize[1] <- sqrt(sum(apply(apply(Y[, 1:2], 2, range), 2, diff)**2))
	lRate <- 2.0 *(eSize[1] +1.0 /eSize[1]) *log(zSize *nnSize)
	clusterExport(cl, c('lRate', 'alpha'), envir = environment())

	# +++ reset number of iterations on workers
	clusterExport(cl, c('epochs', 'iters'), envir = environment())

	# +++ just for testing
	clusterEvalQ(cl, epoch <<- 0)
	
	# start
	t0 <- Sys.time()

	# report starting information
	avgCost <- mean(eCost[ , 1])
	nulL <- ptsne.info(threads, zSize, nnSize, epochs, iters, 0, avgCost, eSize[1], t0, theta)

	for (e in seq(epochs)) {
		# epoch starting time
		te <- Sys.time()
		# resample dataset. Att!! C++ indexes
		I <- sample(seq(bdm$nX)) -1
		# get I&Y chunks
		zChnks(Z.list, Y, I, zBrks)
		# perform partial ptSNE and pool workers epoch-cost
		# NULL returned by thread.rank == 0 is removed by unlist
		eCost[, (e +1)] <- unlist(clusterApply(cl, Z.list, mpi.ztsne))
		# pool partial mappings from workers
		zMap.list <- clusterEvalQ(cl, if (thread.rank != 0) w$zY)
		# restructure global mapping
		updateY(Y, I, zMap.list, zBrks)
		# embedding size
		eSize[(e +1)] <- sqrt(sum(apply(apply(Y[, 1:2], 2, range), 2, diff)**2))
		# learning Rate
		lRate <- 2.0 *(eSize[(e +1)] +1.0 /eSize[(e +1)]) *log(zSize *nnSize)
		# momentum
		if (qDecay)
			alpha <- momentum *(1 -e /scheme$epochs)**2
		else
			alpha <- momentum *(1 -e /scheme$epochs)
		clusterExport(cl, c('lRate', 'alpha'), envir = environment())
		# report status
		if (progress >= 0) {
			avgCost <- mean(eCost[, e+1])
			epoch.info(e, epochs, avgCost, eSize[e +1], t0, te, lRate)
		}
		# security break control
		if (mean(eCost[, e+1]) < 0) break
	}

	bdm$ptsne$Y <- Y
	bdm$ptsne$cost <- as.matrix(eCost[, 1:(epochs +1)])
	bdm$ptsne$size <- eSize

	# report status
	avgCost <- mean(eCost[, e +1])
	bdm$t[['epoch']] <- ptsne.info(threads, zSize, nnSize, epochs, iters, e, avgCost, eSize[e +1], t0, theta)
	#
	return(bdm)
}

# -----------------------------------------------------------------------------
# +++ ptSNE MPI worker functions
# -----------------------------------------------------------------------------

mpi.ztsne <- function(zChnk)
{
	epoch <<- epoch +1	# Att.!! We are using global variable epoch
	if (thread.rank != 0) {
		w$zI <-  zChnk[, 1]
		w$zY <- zChnk[, 2:3]
		mpi_zTSNE(thread.rank, epoch, Xbm@address, Bbm@address, w$zY, w$zI, iters, nnSize, theta, lRate, alpha, gain, is.distance, is.sparse)
	}
}

# -----------------------------------------------------------------------------
# +++ auxiliary functions for ptSNE
# -----------------------------------------------------------------------------

# +++ initial embedding: unity disk
ptsne.init <- function(n, layers)
{
	tht <- runif(n) *2 *pi
	rad <- sqrt(runif(n)) *sqrt(0.5) /2
	Y <- cbind(rad *cos(tht), rad *sin(tht))
	if (layers > 1) {
		for (l in seq(layers -1)) {
			tht <- runif(n) *2 *pi
			rad <- sqrt(runif(n)) *sqrt(0.5) /2
			Y <- cbind(Y, rad *cos(tht), rad *sin(tht))
		}
	}
	return(Y)
}

time.format <- function(time.secs)
{
	time.secs <- round(time.secs, 0)
	if (time.secs > 86400)
	{
		dd <- time.secs %/% 86400
		hh <- (time.secs - dd*86400) %/% 3600
		mm <- (time.secs - dd*86400 - hh*3600) %/% 60
		ft <- paste(formatC(dd, width=2, flag='0'), 'D', formatC(hh, width=2, flag='0'), ':', formatC(mm, width=2, flag='0'), sep='')
	} else
	{
		hh <- time.secs %/% 3600
		mm <- (time.secs - hh*3600) %/% 60
		ss <- (time.secs - hh*3600 - mm*60)
		ft <- paste(formatC(hh, width=2, flag='0'), ':', formatC(mm, width=2, flag='0'), ':', formatC(ss, width=2, flag='0'), sep='')
	}
	return(ft)
}

ptsne.info <- function(threads, zSize, nnSize, epochs, iters, e, avgCost, avgSize, t0, theta)
{
	cat('...')
	cat(' threads ', threads, sep='')
	cat(', size ', zSize, sep='')
	cat(', nnSize ', nnSize, sep='')
	cat(', theta ', theta, sep='')
	cat(', epochs ', epochs, sep='')
	cat(', iters ', iters, sep='')
	cat('\n')

	if (e == 0) {
		cat('--- epoch     ')
		cat('   <Cost>   ')
		cat('   <Size>   ')
		cat('   epSecs ')
		cat(' time2End ')
		cat(' eTime ')
		cat('\n')
	}

	cat('+++ ', formatC(e, width=4, flag='0'), '/', formatC(epochs, width=4, flag='0'), sep='')
	cat(formatC(avgCost, format='e', digits=4, width=12))
	cat(formatC(avgSize, format='e', digits=4, width=12))
	epochSecs <- 0
	if (e > 0){
		runTime <- as.numeric(difftime(Sys.time(), t0, units='secs'))
		epochSecs <- runTime/epochs
		cat('  ', formatC(epochSecs, format='f', digits=4, width=8), sep='')
		cat('  ', time.format(runTime), sep='')

	}
	cat('\n')
	return(epochSecs)
}

epoch.info <- function(e, epochs, avgCost, avgSize, t0, te, lRate)
{
	cat('+++ ', formatC(e, width=4, flag='0'), '/', formatC(epochs, width=4, flag='0'), sep='')
	cat(formatC(avgCost, format='e', digits=4, width=12))
	cat(formatC(avgSize, format='e', digits=4, width=12))
	epchSecs <- as.numeric(difftime(Sys.time(), te, units='secs'))
	cat('  ', formatC(epchSecs, format='f', digits=4, width=8), sep='')
	runTime <- as.numeric(difftime(Sys.time(), t0, units='secs'))
	tm2End <- runTime/e * (epochs-e)
	cat('  ', time.format(tm2End), sep='')
	cat('  ', format(Sys.time()+tm2End, '%H:%M'), sep='')
	cat('  ', lRate)
	cat('\n')
}
