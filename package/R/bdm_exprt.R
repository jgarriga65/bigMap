# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# +++ Export input-data (normalize = F for bdm.pakde())
# -----------------------------------------------------------------------------

get.Xbm.desc <- function()
{
	if (thread.rank == thread.hldr) {
		if (normalize) centerScale(Xbm@address, is.distance, is.sparse)
		describe(Xbm)
	}
}

Xdata.exp <- function(cl, inp.data, is.distance, is.sparse, normalize)
{
	clusterExport(cl, c('is.distance', 'is.sparse', 'normalize'), envir = environment())
	if (attr(cl[[1]], 'class') == 'SOCKnode')
	{
		if (class(inp.data) == 'character') {
			cat('+++ loading ', inp.data, ' \n', sep = '')
			dataFile <- inp.data
			inp.data <- as.matrix(read.csv(dataFile))
		}
		else if (class(inp.data) == 'big.matrix') {
			inp.data <- inp.data[, ]
		}
		Xbm <- as.big.matrix(t(inp.data), type = 'double')
		# by default normalize input-data (avoids unclear numerical problems)
		if (normalize) centerScale(Xbm@address, is.distance, is.sparse)
		Xbm.dsc <- describe(Xbm)
		# export big matrix descriptor to workers
		clusterExport(cl, c('Xbm.dsc'), envir = environment())
		# attach big matrix to workers
		clusterEvalQ(cl, Xbm <- attach.big.matrix(Xbm.dsc))
	}
	else if (class(inp.data) == 'character')
	{
		# attach big.matrix data to holders
		cat('+++ loading ', inp.data, ' \n', sep = '')
		dataFile <- inp.data
		clusterExport(cl, c('dataFile', 'is.distance'), envir = environment())
		cl.Xdsc <- clusterEvalQ(cl,
			if (thread.rank == thread.hldr) {
				Xbm <- as.big.matrix(t(as.matrix(read.csv(dataFile))), type = 'double')
			}
		)
		# get big.matrix backing.file descriptors from holders
		# by default normalize input-data (avoids unclear numerical problems)
		cl.Xdsc <- clusterCall(cl, get.Xbm.desc)
		# export shared-memory descriptors
		clusterExport(cl, c('cl.Xdsc'), envir = environment())
		# attach big matrix to workers
		nulL <- clusterEvalQ(cl,
			if (thread.rank != thread.hldr){
				Xbm <- attach.big.matrix(cl.Xdsc[[thread.hldr]])
			}
		)
	}
	else
	{
		if (class(inp.data) == 'big.matrix') inp.data <- inp.data[, ]
		f <- tName.get('X')
		Xbf <- as.big.matrix(inp.data, type='double', backingpath = f$path, backingfile = f$bin, descriptorfile = f$desc)
		Xbf.dsc <- describe(Xbf)
		clusterExport(cl, c('Xbf.dsc'), envir = environment())
		# attach big.matrix backing.file to holders
		nulL <- clusterEvalQ(cl,
			if (thread.rank == thread.hldr) {
				Xhl <- attach.big.matrix(Xbf.dsc)
				# Att!! Xhl is a backed file.
				# Workers attached to Xhl read from a backed file.
				# This is extremely slow !!! Must load it to shared memory.
				Xbm <- as.big.matrix(t(Xhl[ , ]), type='double')
				rm(Xhl)
			})
		# get big.matrix backing.file descriptors from holders
		# by default normalize input-data (avoids unclear numerical problems)
		cl.Xdsc <- clusterCall(cl, get.Xbm.desc)
		# export shared-memory descriptors
		clusterExport(cl, c('cl.Xdsc'), envir = environment())
		# attach big matrix to workers
		nulL <- clusterEvalQ(cl,
			if (thread.rank != thread.hldr){
				Xbm <- attach.big.matrix(cl.Xdsc[[thread.hldr]])
			})
		# remove backing file
		unlink(paste(f$path, f$bin, sep = '/'))
		unlink(paste(f$path, f$desc, sep = '/'))
	}
}

Xdata.dim <- function(cl)
{
	# Att.!! X transposed
	nX <- unlist(clusterEvalQ(cl, if (thread.rank == 1) ncol(Xbm)))
	mX <- unlist(clusterEvalQ(cl, if (thread.rank == 1) nrow(Xbm)))
	return(c(nX, mX))
}

# -----------------------------------------------------------------------------
# +++ Export beta to workers
# -----------------------------------------------------------------------------

Xbeta.exp <- function(cl, B)
{
	if (attr(cl[[1]], 'class') == 'SOCKnode')
	{
		# define Bbm big.matrix (betas)
		Bbm <- as.big.matrix(t(B), type = 'double')
		Bbm.dsc <- describe(Bbm)
		# export big matrix descriptor to workers
		clusterExport(cl, c('Bbm.dsc'), envir = environment())
		# attach big matrix to workers
		clusterEvalQ(cl, Bbm <- attach.big.matrix(Bbm.dsc))
	}
	else
	{
		f <- tName.get('B')
		Bbf <- as.big.matrix(t(B), type='double', backingpath = f$path, backingfile = f$bin, descriptorfile = f$desc)
		Bbf.dsc <- describe(Bbf)
		clusterExport(cl, c('Bbf.dsc'), envir = environment())
		# attach big.matrix backing.file to holders
		# get big.matrix backing.file descriptors from holders
		cl.Bdsc <- clusterEvalQ(cl,
			if (thread.rank == thread.hldr) {
				Bhl <- attach.big.matrix(Bbf.dsc)
				Bbm <- as.big.matrix(Bhl[, ], type='double')
				rm(Bhl)
				describe(Bbm)
			}
		)
		# export shared memory descriptors
		clusterExport(cl, c('cl.Bdsc'), envir = environment())
		# attach big matrix to workers
		nulL <- clusterEvalQ(cl,
			if (thread.rank != thread.hldr){
				Bbm <- attach.big.matrix(cl.Bdsc[[thread.hldr]])
			}
		)
		# remove backing file
		unlink(paste(f$path, f$bin, sep = '/'))
		unlink(paste(f$path, f$desc, sep = '/'))
	}
}

# -----------------------------------------------------------------------------
# +++ Export OUTPUT data
# -----------------------------------------------------------------------------

Ydata.exp <- function(cl, Y)
{
	if (attr(cl[[1]], 'class') == 'SOCKnode')
	{
		# define Ybm big.matrix (Y)
		Ybm <- as.big.matrix(Y, type = 'double')
		Ybm.dsc <- describe(Ybm)
		# export big matrix descriptor to workers
		clusterExport(cl, c('Ybm.dsc'), envir = environment())
		# attach big matrix to workers
		clusterEvalQ(cl, Ybm <- attach.big.matrix(Ybm.dsc))
	}
	else
	{
		# attach big.matrix data to holders
		f <- tName.get('Y')
		Ybf <- as.big.matrix(Y, type='double', backingpath = f$path, backingfile = f$bin, descriptorfile = f$desc)
		Ybf.dsc <- describe(Ybf)
		# attach big.matrix backing file to holders
		# get big.matrix descriptors from holders
		clusterExport(cl, c('Ybf.dsc'), envir = environment())
		cl.Ydsc <- clusterEvalQ(cl,
			if (thread.rank == thread.hldr) {
				Yhl <- attach.big.matrix(Ybf.dsc)
				Ybm <- as.big.matrix(Yhl[, ], type='double')
				rm(Yhl)
				describe(Ybm)
			}
		)
		# attach big matrix to workers
		clusterExport(cl, c('cl.Ydsc'), envir = environment())
		nulL <- clusterEvalQ(cl,
			if (thread.rank != thread.hldr){
				Ybm <- attach.big.matrix(cl.Ydsc[[thread.hldr]])
			}
		)
		# remove backing file
		unlink(paste(f$path, f$bin, sep = '/'))
		unlink(paste(f$path, f$desc, sep = '/'))
	}
}
