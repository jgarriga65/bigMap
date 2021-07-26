
# --------------------------------------------------------------------------------
# +++ BHt-SNE and refinement t-SNE with naive parallelization
# --------------------------------------------------------------------------------

bdm.bhtsne <- function(dSet.data, dSet.name = NULL, is.distance = F, is.sparse = F, normalize = F, dSet.labels = NULL, ppx = 100, xppx = 3.0, iters = 100, theta = .5, eSizeX = 2.0 *sqrt(2), exgg = 1, threads = 4)
{
	m <- bdm.init(dSet.data, dSet.name = dSet.name, is.distance = is.distance, is.sparse = is.sparse, normalize = normalize, ppx = ppx, xppx = xppx, dSet.labels = dSet.labels, threads = threads, mpi.cl = NULL)
	m$ppx <- m$ppx[[1]]

	m.list <- bdm.rtsne(dSet.data, m, ppx = ppx, xppx = xppx, iters = iters, theta = theta, eSizeX = eSizeX, exgg = exgg, threads = threads)
	return(m.list[[2]])
}

bdm.rtsne <- function(dSet.data, m, ppx, xppx = 3.0, iters = 100, theta = .5, eSizeX = 2.0, exgg = 1, threads = 4)
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
	if (m$ppx$ppx != ppx || m$ppx$xppx != xppx) {
		m$ppx <- beta.get(cl, ppx, xppx)
	}
	cat('+++ exporting betas \n')
	Xbeta.exp(cl, m$ppx$B)
	# +++ define thread chunks
	# Att!! nnSize refers here to the whole dataset, thus max is (nX -1)
	nX <- m$nX
	lognX <- log(nX * (nX -1))
	nnSize <- min((m$ppx$ppx *m$ppx$xppx), (nX -1))
	# +++ initial embedding
	if (is.null(m$ptsne$Y)){
		Y <- ptsne.init(m$nX, 1)
	} else {
		Y <- m$ptsne$Y[, 1:2]
	}
	mY <- 2
	clusterExport(cl, c('nX', 'mY', 'nnSize'), envir = environment())
	# +++ thread initialization
	clusterCall(cl, thread.init)
	# +++ compute thread affinity matrices
	m$t$affMtx <- system.time({
		cat('+++ computing thread affMtx \n')
		rawP <- unlist(clusterCall(cl, thread.affMtx))
		print(rawP)
		rawP <- sum(rawP)
		cat('... sumP', rawP, '/', nX, formatC(rawP /nX, format = 'e', digits = 4, width = 12), '\n')
		sumP <- rawP /exgg
		clusterExport(cl, c('sumP'), envir = environment())
	})
	print(m$t$affMtx)
	# +++ start bh-tsne
	# m$progress <- list()
	t0 <- Sys.time()
	info.head()
	infoRate <- 1
	itCost <- rep(0, iters)
	itSize <- rep(0, iters)
	# +++ BHt-SNE parameters
	eSize_ <- sqrt(sum(apply(apply(Y, 2, range), 2, diff)**2))
	eSizeX <- eSizeX *sqrt(2)
	lRate <- 2.0 *eSizeX *log2(nX *nnSize)
	alpha <- 0.5
	clusterExport(cl, c('lRate', 'theta', 'alpha'), envir = environment())
	shiftIt <- 25
	m$t$bh <- system.time({
		# +++ iterate
		for (it in seq(iters)) {
			t1 <- system.time({
				# +++ export current embedding
				# (no need to transpose as each thread we'll use a local copy!!!)
				eSize <- sqrt(sum(apply(apply(Y, 2, range), 2, diff)**2))
				lRate <- eSizeX *(eSize /eSize_ +eSize_ /eSize) *log2(nX *nnSize)
				clusterExport(cl, c('lRate'), envir = environment())
				if (it < iters /2) {
					alpha <- 0.5 +0.6 *it /iters
					clusterExport(cl, c('alpha'), envir = environment())
				}
				itSize[it] <- eSize
				Ydata.exp(cl, Y)
				# +++ repulsice forces
				sumQ <- sum(unlist(clusterCall(cl, thread.repFget)))
				clusterExport(cl, c('sumQ'), envir = environment())
				# +++ update embedding
				Y <- matrix(unlist(clusterCall(cl, thread.itrRun)), ncol = 2, byrow = T)
				eCost <- sum(unlist(clusterEvalQ(cl, if (thread.rank != 0) zCost)))
				itCost[it] <- (log(sumQ) -eCost /rawP) /lognX
				nulL <- clusterCall(cl, thread.rmYbm)
				# +++ check regime shift
				if (it == shiftIt) {
					infoRate <- 25
					# if (abs(itCost[it] -itCost[(it -10)]) /itCost[(it -10)] < 1e-4) shiftIt <- it
				}
			})
			# m$progress[[it]] <- list(epoch = it, Y = Y)
			if (it %% infoRate == 0) iter.info(it, iters, t0, t1, itCost[it], itSize[it], lRate)
			# if (it == shiftIt) {
			# 	sumP <- sumP *exgg
			# 	clusterExport(cl, c('sumP'), envir = environment())
			# }
		}
	})
	if (it %%infoRate != 0) {
		cat('... \n')
		iter.info(it, iters, t0, t1, eCost, eSize, lRate)
	}
	nulL <- clusterCall(cl, thread.rmAll)
	print(m$t$bh)
	# +++ cluster stop
	cluster.stop(cl)
	m$ptsne <- list(threads = 1, layers = 1, bh.threads = threads, iters = iters, theta = theta, alpha = alpha, lRate = lRate, exgg = exgg, shitIt = shiftIt)
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
		Ubm <<- as.big.matrix(matrix(rep(0, (z.end -z.ini) *mY), nrow = nnSize), type = 'double')
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

iter.info <- function(it, iters, t0, t1, eCost, eSize, lRate)
{
	cat('+++ ', formatC(it, width=3, flag='0'), '/', formatC(iters, width=3, flag='0'), sep='')
	cat(formatC(eCost, format='e', digits=4, width=12))
	cat(formatC(eSize, format='e', digits=4, width=12))
	cat('  ', formatC(t1[3], format='f', digits=4, width=8), sep='')
	runTime <- as.numeric(difftime(Sys.time(), t0, units='secs'))
	tm2End <- runTime /it * (iters -it)
	cat('   ', time.format(tm2End), sep='')
	cat('  ', format(Sys.time()+tm2End, '%H:%M'), sep='')
	cat('  ', lRate)
	cat('\n')
}
