% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bdm_main.R
\name{bdm.ptsne}
\alias{bdm.ptsne}
\title{Parallelized t-SNE (ptSNE)}
\usage{
bdm.ptsne(
  data,
  bdm,
  theta = 0.5,
  Y.init = NULL,
  mpi.cl = NULL,
  threads = 4,
  layers = 2,
  info = 0
)
}
\arguments{
\item{data}{Input data (a matrix, a big.matrix or a .csv file name).}

\item{bdm}{A \var{bdm} data mapping instance.}

\item{theta}{Accuracy/speed trade-off factor, a value between 0.33 and 0.8. (Default value is \var{theta = 0.0}). If \var{theta < 0.33} the algorithm uses the exact computation of the gradient. The closer is this value to 1 the faster is the computation but the coarser is the approximation of the gradient.}

\item{Y.init}{A \var{n *2 *layers} matrix with initial mapping positions. (By default \var{Y.init=NULL} will use random initial positions).}

\item{mpi.cl}{MPI (inter-node parallelization) cluster as generated by \var{bdm.mpi.start()}. (By default \var{mpi.cl = NULL} a 'SOCK' (intra-node parallelization) cluster is generated).}

\item{threads}{Number of parallel threads (according to data size and hardware resources, \var{i.e.} number of cores and available memory. Default value is \var{threads = 4}).}

\item{layers}{Number of layers (\var{minimum} 2, \var{maximum} the number of threads). Default value is \var{layers = 2}.}

\item{info}{Output information: 1 yields inter-round results, 0 disables intermediate results. Default value is \var{info = 0}.}
}
\value{
A \var{bdm} data mapping instance.
}
\description{
Starts the parallelized t-SNE algorithm (pt-SNE). This is the first step of the mapping protocol.
}
\examples{

# --- load example dataset
bdm.example()
# --- perform ptSNE
\dontrun{
# --- run ptSNE
m <- bdm.ptsne(ex$data, ex$map, threads = 10, layers = 2)
# --- plot the Cost function
bdm.cost(m)
# --- plot ptSNE output
bdm.ptsne.plot(m, class.lbls = ex$labels)
}
}
