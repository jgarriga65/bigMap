% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bdm_mtsne.R
\name{bdm.mtsne}
\alias{bdm.mtsne}
\title{Multi-core t-SNE (mtSNE)}
\usage{
bdm.mtsne(
  dSet.data,
  is.distance = F,
  is.sparse = F,
  ppx = 100,
  theta = 0.5,
  iters = 250,
  mpi.cl = NULL,
  threads = 4,
  infoRate = 25
)
}
\arguments{
\item{dSet.data}{A \var{data.frame} or \var{matrix} with raw input-data. The dataset must not have duplicated rows.}

\item{is.distance}{Default value is \var{is.distance = FALSE}. TRUE indicates that raw data is a distance matrix.}

\item{is.sparse}{Default value is \var{is.sparse = FALSE}. TRUE indicates that the raw data is a sparse matrix.}

\item{ppx}{The value of perplexity to compute similarities.}

\item{theta}{Accuracy/speed trade-off factor, a value between 0.33 and 0.8. Default value is \var{theta = 0.5}. If \var{theta < 0.33} the algorithm uses the exact computation of the gradient. The closer it is this value to 1 the faster the computation and the coarser the approximation of the gradient.}

\item{iters}{Number of iters/epoch. Default value is \var{iters = 250}.}

\item{mpi.cl}{An MPI (inter-node parallelization) cluster as generated by \var{bdm.mpi.start()}. By default \var{mpi.cl = NULL}, i.e. a 'SOCK' (intra-node parallelization) cluster is generated.}

\item{threads}{Number of parallel threads (according to data size and hardware resources, i.e. number of cores and available memory). Default value is \var{threads = 4}.}

\item{infoRate}{Number of epochs to show output information. Default value is \var{infoRate = 25}.}
}
\value{
A \var{bdm} data mapping instance.
}
\description{
Starts the multi-core t-SNE (mtSNE) algorithm.
}
\examples{

# --- load example dataset
bdm.example()
\dontrun{
# --- perform mtSNE
m <- bdm.mtsne(ex$data, ex$map, ppx = 250, iters = 250, threads = 4)
# --- plot the Cost function
bdm.cost(m)
# --- plot mtSNE output (use bdm.ptsne.plot() function)
bdm.ptsne.plot(m)
}
}
