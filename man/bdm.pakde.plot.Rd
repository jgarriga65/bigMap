% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bdm_plot.R
\name{bdm.pakde.plot}
\alias{bdm.pakde.plot}
\title{Plot paKDE (density landscape)}
\usage{
bdm.pakde.plot(bdm, pakde.pltt = NULL, pakde.lvls = 16, layer = 1)
}
\arguments{
\item{bdm}{A \var{bdm} data mapping instance, or a list of them to make a comparative plot.}

\item{pakde.pltt}{A colour palette to show levels in the paKDE plot. By default (\var{pakde.pltt = NULL}) the default palette is used.}

\item{pakde.lvls}{The number of levels of the density heat-map (16 by default).}

\item{layer}{The ptSNE output layer. Default value is \var{layer = 1}.}
}
\value{
None.
}
\description{
Plot paKDE (density landscape)
}
\examples{

bdm.example()
bdm.pakde.plot(ex$map)
}
