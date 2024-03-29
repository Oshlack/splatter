% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kersplat-simulate.R
\name{kersplatGenNetwork}
\alias{kersplatGenNetwork}
\title{Generate Kersplat gene network}
\usage{
kersplatGenNetwork(params, verbose)
}
\arguments{
\item{params}{KersplatParams object containing simulation parameters.}

\item{verbose}{logical. Whether to print progress messages}
}
\value{
KersplatParams object with gene network
}
\description{
Generate a gene network for the Kersplat simulation
}
\details{
Currently a very simple approach is used which needs to be improved. A
network is generated using the \code{\link[igraph]{sample_forestfire}}
function and edge weights are sampled from a standard normal distribution.
}
\keyword{internal}
