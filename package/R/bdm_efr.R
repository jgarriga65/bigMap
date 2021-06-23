
# -----------------------------------------------------------------------------
# +++ Export matrix of P,Q normalization factors
# -----------------------------------------------------------------------------

zQ.exp <- function(cl, threads)
{
	if (attr(cl[[1]], 'class') == 'SOCKnode')
	{
		bigZ <- as.big.matrix(matrix(rep(0, threads), ncol = 1), type = 'double')
		bigZ.dsc <- describe(bigZ)
		# export big matrix descriptor to workers
		clusterExport(cl, c('bigZ.dsc'), envir = environment())
		# attach big matrix to workers
		clusterEvalQ(cl, bigZ <- attach.big.matrix(bigZ.dsc))
	}
	else
	{
		# attach big.matrix data to holders
		f <- tName.get('Z')
		bckZ <- as.big.matrix(matrix(rep(0, threads), ncol = 1), type='double', backingpath = f$path, backingfile = f$bin, descriptorfile = f$desc)
		bckZ.dsc <- describe(bckZ)
		# attach big.matrix backing file to holders
		# get big.matrix descriptors from holders
		clusterExport(cl, c('bckZ.dsc'), envir = environment())
		cl.Zdsc <- clusterEvalQ(cl,
			if (thread.rank == thread.hldr) {
				hldZ <- attach.big.matrix(bckZ.dsc)
				bigZ <- as.big.matrix(as.matrix(hldZ[ , ]), type='double')
				rm(hldZ)
				describe(bigZ)
			}
		)
		# attach big matrix to workers
		clusterExport(cl, c('cl.Zdsc'), envir = environment())
		nulL <- clusterEvalQ(cl,
			if (thread.rank != thread.hldr){
				bigZ <- attach.big.matrix(cl.Zdsc[[thread.hldr]])
			}
		)
		# remove backing file
		unlink(paste(f$path, f$bin, sep = '/'))
		unlink(paste(f$path, f$desc, sep = '/'))
	}
}

# --------------------------------------------------------------------------------
# +++ dataset		n		ppx		theta	cl		np	GB	node	time
# ... technosperm 	63931			0.0		MPI		10	4	101		3.80 1.90 3.30
# ... technosperm 	63931			0.0		SCKT	10	4	105		0.05 0.25 1.60
# ... macosko15 	44808			0.0		MPI		20	4	100		9.26 5.36 3.29
# ... macosko15 	44808			0.0		SCKT	10	4	111		0.05 1.12 1.50
# ... P1G			1000000	500		0.33	SCKT	15	8	111		0.13 2.68 6.15
# ... P1G			1000000	50		0.33	SCKT	15	8	111		0.09 2.68 3.25 
# --------------------------------------------------------------------------------

bdm.efr <- function(dSet.data, m.list, ppx = ppx, iters = 100, lRate = NULL, theta = .0, exgg = 1, threads = 4, mpi.cl = NULL)
{
	m <- m.list[[1]]
	# +++ start cluster
	cl <- cluster.start(threads, mpi.cl)
	# +++ export input data (if using mpi.cl it might have been already exported)
	if (is.null(mpi.cl) || !is.null(dSet.data)) {
		Xdata.exp(cl, dSet.data, m$is.distance, m$is.sparse, m$normalize)
	}
	# +++ compute/export betas
	if (class(ppx) == 'list') {
		m$ppx <- ppx
	}
	else if (ppx != m$ppx$ppx) {
		m$ppx <- beta.get(cl, ppx, m$ppx$xppx)
	}
	Xbeta.exp(cl, m$ppx$B)
	# +++ define thread chunks
	nX <- m$nX
	# Att!! nnSize refers here to the whole dataset
	nnSize <- min((m$ppx$ppx *m$ppx$xppx), (nX -1))
	clusterExport(cl, c('nX', 'nnSize'), envir = environment())
	clusterCall(cl, thread.chunk)
	# +++ compute thread affinity matrices
	cat('+++ computing thread affMtx \n')
	t <- system.time({
		zP <- unlist(clusterCall(cl, thread.affMtx))
		print(zP /exgg)
		zP <- sum(zP)/exgg
		clusterExport(cl, c('zP'), envir = environment())
	})
	print(t)
	# +++ ptSNE parameters
	if (is.null(lRate)) lRate <- nX /16
	clusterExport(cl, c('lRate', 'theta'), envir = environment())
	# +++ embedding final compression
	# m$progress <- list()
	m$t$efr <- system.time({
		# +++ output affinities (Q) normalization
		# cat('+++ export Q normalization vector \n', sep='')
		# zQ.exp(cl, threads)
		# +++ initial embedding (layer = 1)
		Y <- m$ptsne$Y[, 1:2]
		# +++ iterate
		for (it in seq(iters)) {
			# +++ export embedding (no need to transpose as each thread we'll use a local copy!!!)
			t1 <- system.time(Ydata.exp(cl, Y))
			t2 <- system.time({
				zQ <- unlist(clusterCall(cl, thread.qNorm))
				zQ <- sum(zQ)
				clusterExport(cl, c('zQ'), envir = environment())
			})
			t3 <- system.time(Y <- matrix(unlist(clusterCall(cl, thread.itrRun)), ncol = 2, byrow = T))
			nulL <- clusterCall(cl, thread.rmYbm)
			# m$progress[[it]] <- list(epoch = it, Y = Y)
			cat('+++ iteration ', it, zQ, t1[3], t2[3], t3[3], '\n')
		}
	})
	print(m$t$efr)
	# +++ cluster stop
	if (is.null(mpi.cl)) cluster.stop(cl)
	m$ptsne$Y <- Y
	m$ptsne$efrIt <- iters
	m.list[[(length(m.list) +1)]] <- m
	return(m.list)
}

thread.chunk <- function()
{
	if (thread.rank != 0) {
		breaks <- c(sapply(1:threads, function(rank) (rank -1) *(nX +1.0) %/%threads), nX)
		# C++ indexes
		z.ini <<- breaks[thread.rank];
		z.end <<- breaks[thread.rank +1];
	}
}

thread.affMtx <- function()
{
	if (thread.rank != 0) {
		Pbm <<- as.big.matrix(matrix(rep(0, (z.end -z.ini) *nnSize), nrow = nnSize), type = 'double')
		Wbm <<- as.big.matrix(matrix(rep(0, (z.end -z.ini) *nnSize), nrow = nnSize), type = 'integer')
		zP <- thread_affMtx(z.ini, z.end, Xbm@address, is.distance, is.sparse, Bbm@address, Pbm@address, Wbm@address)
		zP
	}
}

thread.qNorm <- function()
{
	if (thread.rank != 0) {
		zQ <- thread_qNorm(z.ini, z.end, Ybm@address, theta)
		zQ
	}
}

thread.itrRun <- function()
{
	if (thread.rank != 0) {
		thread_itrRun(z.ini, z.end, Pbm@address, Wbm@address, Ybm@address, zP, zQ, lRate, theta)
	}
}

thread.rmYbm <- function()
{
	if (thread.rank != 0) rm(Ybm)
}
