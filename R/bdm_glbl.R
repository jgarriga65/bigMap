# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# cluster configuration
varNames <- c('threads', 'thread.rank', 'thread.hldr')

# shared data
varNames <- c(varNames, 'Xbm', 'Xbf.dsc', 'cl.Xdsc')

# shared Betas
varNames <- c(varNames, 'Bbm', 'Bbf.dsc', 'cl.Bdsc')

# shared embedding
varNames <- c(varNames, 'zY', 'Ybm')

# shared row indexes
varNames <- c(varNames, 'Ibm')

# computation of betas
varNames <- c(varNames, 'itr', 'tol', 'ppx')

# thread affinities
varNames <- c(varNames, 'z.ini', 'z.end')

# ptSNE parameters
varNames <- c(varNames, 'layers', 'iters', 'theta', 'alpha', 'lRate', 'is.distance', 'is.sparse', 'info')

# thread ptSNE variables
varNames <- c(varNames, 'sumP', 'sumQ', 'epoch', 'zSize', 'nnSize')

# mtSNE parameters
varNames <- c(varNames, 'nX', 'mY', 'Pbm', 'Wbm', 'Rbm', 'Ubm', 'Gbm')

# pakde function
varNames <- c(varNames, 'thread.pakde', 'pakde.grid')

# dMap function
varNames <- c(varNames, 'pakde', 'Y', 'L')

# kNP + hlCorr
varNames <- c(varNames, 'k.max', 'sampling', 'zSampleSize') 

# Main settings (do not change !!)

xppx <- 3
normalize <- T

momentum <- 0.8
gain <- 2.0
qDecay <- T

# Global variables declaration

if (getRversion() >= '2.15.1')
    utils::globalVariables(varNames, 'bigMap', add = TRUE)
