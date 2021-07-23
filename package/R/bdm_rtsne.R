
# --------------------------------------------------------------------------------
# +++ BHt-SNE and refinement t-SNE with naive parallelization
# --------------------------------------------------------------------------------

bdm.bhtsne <- function(dSet.data, dSet.name = NULL, is.distance = F, is.sparse = F, normalize = F, dSet.labels = NULL, ppx = 100, iters = 100, theta = .0, lRate = NULL, exgg = 1, threads = 4)
{
	m <- bdm.init(dSet.data, dSet.name = dSet.name, is.distance = is.distance, is.sparse = is.sparse, normalize = normalize, ppx = ppx, xppx = 3.0, dSet.labels = dSet.labels, threads = threads, mpi.cl = NULL)
	m$ppx <- m$ppx[[1]]

	Y.init <- ptsne.init(m$nX, 1)
	m$ptsne <- list(threads = 1, layers = 1, bhthreads = threads, iters = iters, theta = theta, lRate = lRate, exgg = exgg, shiftIt = iters /4, Y = Y.init)

	m.list <- bdm.rtsne(dSet.data, m, ppx = ppx, iters = iters, theta = theta, lRate = lRate, exgg = exgg, threads = threads)
	return(m.list[[2]])
}

bdm.rtsne <- function(dSet.data, m, ppx = ppx, iters = 100, theta = .5, lRate = NULL, exgg = 1, threads = 4)
{
	m.list <- list(m)
	# +++ start cluster
	cl <- cluster.start(threads, NULL)
	# +++ export input data (if using mpi.cl it might have been already exported)
	cat('+++ exporting data \n')
	if (!is.null(dSet.data)) {
		Xdata.exp(cl, dSet.data, m$is.distance, m$is.sparse, m$normalize)
	}
	# +++ compute/export betas
	# Att!! xppx must be 3 here !!!
	if (m$ppx$ppx != ppx || m$ppx$xppx != 3) {
		m$ppx <- beta.get(cl, ppx, 3)
	}
	cat('+++ exporting betas \n')
	Xbeta.exp(cl, m$ppx$B)
	# +++ define thread chunks
	# Att!! nnSize refers here to the whole dataset, thus max is (nX -1)
	nX <- m$nX
	lognX <- log(nX * (nX -1))
	nnSize <- min((m$ppx$ppx *m$ppx$xppx), (nX -1))
	# +++ initial embedding (layer = 1)
	Y <- m$ptsne$Y[, 1:2]
	mY <- 2
	clusterExport(cl, c('nX', 'mY', 'nnSize'), envir = environment())
	# +++ thread initialization
	clusterCall(cl, thread.init)
	# +++ compute thread affinity matrices
	cat('+++ computing thread affMtx \n')
	t <- system.time({
		rawP <- unlist(clusterCall(cl, thread.affMtx))
		print(rawP /exgg)
		rawP <- sum(rawP)
		sumP <- rawP /exgg
		clusterExport(cl, c('sumP'), envir = environment())
	})
	print(t)
	# +++ ptSNE parameters
	eRange <- sqrt(sum(apply(apply(Y, 2, range), 2, diff)**2))
	if (is.null(lRate)) {
		# lRate <- (eRange /2) +(2 /eRange) *log2(nX *nnSize)
		lRate <- (nX -1) /exgg
	}
	alpha <- 0.5
	clusterExport(cl, c('lRate', 'theta', 'alpha'), envir = environment())
	# +++ start bh-tsne
	# m$progress <- list()
	t0 <- Sys.time()
	info.head()
	infoRate <- 1
	itCost <- rep(0, iters)
	itSize <- rep(0, iters)
	shiftIt <- iters /4
	m$t$bh <- system.time({
		# +++ iterate
		for (it in seq(iters)) {
			t1 <- system.time({
				# +++ export current embedding
				# (no need to transpose as each thread we'll use a local copy!!!)
				itSize[it] <- sqrt(sum(apply(apply(Y, 2, range), 2, diff)**2))
				# lRate <- itSize[it] /2 *log2(nX *nnSize)
				# lRate <- (itSize[it] /2) +(2 /itSize[it]) *log2(nX *nnSize)
				# clusterExport(cl, c('lRate'), envir = environment())
				Ydata.exp(cl, Y)
				# +++ repulsice forces
				sumQ <- sum(unlist(clusterCall(cl, thread.repFget)))
				clusterExport(cl, c('sumQ'), envir = environment())
				# +++ update embedding
				Y <- matrix(unlist(clusterCall(cl, thread.itrRun)), ncol = 2, byrow = T)
				Cost <- sum(unlist(clusterEvalQ(cl, if (thread.rank != 0) zCost)))
				itCost[it] <- (log(sumQ) -Cost /rawP) /lognX
				nulL <- clusterCall(cl, thread.rmYbm)
				# +++ check regime shift
				if (it > 10) {
					if (abs(itCost[it] -itCost[(it -10)]) /itCost[(it -10)] < 1e-4) shiftIt <- it
				}
			})
			# m$progress[[it]] <- list(epoch = it, Y = Y)
			if (it %% infoRate == 0) iter.info(it, iters, t0, t1, itCost[it], itSize[it], sumQ)
			if (it == shiftIt) {
				alpha <- 0.8
				sumP <- sumP *exgg
				clusterExport(cl, c('alpha', 'sumP'), envir = environment())
				infoRate <- 10
			}
		}
	})
	if (it %%infoRate != 0) iter.info(it, iters, t0, t1, Cost, eRange)
	nulL <- clusterCall(cl, thread.rmAll)
	print(m$t$bh)
	# +++ cluster stop
	cluster.stop(cl)
	m$ptsne$threads <- 1
	m$ptsne$layers <- 1
	m$ptsne$bhthreads <- threads
	m$ptsne$iters <- iters
	m$ptsne$lRate <- lRate
	m$ptsne$exgg <- exgg
	m$ptsne$shiftIt <- shiftIt
	m$ptsne$Y <- Y
	m$ptsne$cost <- matrix(itCost, ncol = 1)
	m$ptsne$size <- itSize
	m.list[[(length(m.list) +1)]] <- m
	return(m.list)
}

thread.init <- function()
{
	if (thread.rank != 0) {
		breaks <- c(sapply(1:threads, function(rank) (rank -1) *(nX +1.0) %/%threads), nX)
		# C++ indexes
		z.ini <<- breaks[thread.rank];
		z.end <<- breaks[thread.rank +1];
		# local matrices
		Pbm <<- as.big.matrix(matrix(rep(0, (z.end -z.ini) *nnSize), nrow = nnSize), type = 'double')
		Wbm <<- as.big.matrix(matrix(rep(0, (z.end -z.ini) *nnSize), nrow = nnSize), type = 'integer')
		Rbm <<- as.big.matrix(matrix(rep(0, (z.end -z.ini) *mY), nrow = nnSize), type = 'double')
		Ubm <<- as.big.matrix(matrix(rep(1, (z.end -z.ini) *mY), nrow = nnSize), type = 'double')
		Gbm <<- as.big.matrix(matrix(rep(1, (z.end -z.ini) *mY), nrow = nnSize), type = 'double')
		# thread cost
		zCost <<- .0
	}
}

thread.affMtx <- function()
{
	if (thread.rank != 0) {
		thread_affMtx(z.ini, z.end, Xbm@address, is.distance, is.sparse, Bbm@address, Pbm@address, Wbm@address)
	}
}

thread.repFget <- function()
{
	if (thread.rank != 0) {
		thread_repF(z.ini, z.end, Ybm@address, theta, Rbm@address)
	}
}

thread.itrRun <- function()
{
	if (thread.rank != 0) {
		zRet <- thread_iter(z.ini, z.end, Pbm@address, Wbm@address, Ybm@address, sumP, sumQ, Rbm@address, Ubm@address, Gbm@address, lRate, alpha)
		zCost <<- zRet$zCost
		zRet$newY
	}
}

thread.rmYbm <- function()
{
	if (thread.rank != 0) rm(Ybm)
}

thread.rmAll <- function()
{
	if (thread.rank != 0) {
		rm(Pbm)
		rm(Wbm)
		rm(Rbm)
		rm(Ubm)
	}
}

info.head <- function()
{
	cat('--- iter     ')
	cat('   <Qlty>   ')
	cat('   <Size>   ')
	cat(' <itSecs> ')
	cat(' time2End ')
	cat(' eTime ')
	cat('\n')
}

iter.info <- function(it, iters, t0, t1, Cost, eRange, sumQ)
{
	cat('+++ ', formatC(it, width=3, flag='0'), '/', formatC(iters, width=3, flag='0'), sep='')
	cat(formatC(Cost, format='e', digits=4, width=12))
	cat(formatC(eRange, format='e', digits=4, width=12))
	cat('  ', formatC(t1[3], format='f', digits=4, width=8), sep='')
	runTime <- as.numeric(difftime(Sys.time(), t0, units='secs'))
	tm2End <- runTime /it * (iters -it)
	cat('   ', time.format(tm2End), sep='')
	cat('  ', format(Sys.time()+tm2End, '%H:%M'), sep='')
	cat('  ', sumQ)
	cat('\n')
}
