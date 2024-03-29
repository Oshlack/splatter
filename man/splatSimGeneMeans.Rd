% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splat-simulate.R
\name{splatSimGeneMeans}
\alias{splatSimGeneMeans}
\title{Simulate gene means}
\usage{
splatSimGeneMeans(sim, params)
}
\arguments{
\item{sim}{SingleCellExperiment to add gene means to.}

\item{params}{SplatParams object with simulation parameters.}
}
\value{
SingleCellExperiment with simulated gene means.
}
\description{
Simulate gene means from a gamma distribution. Also simulates outlier
expression factors. Genes with an outlier factor not equal to 1 are replaced
with the median mean expression multiplied by the outlier factor.
}
\keyword{internal}
