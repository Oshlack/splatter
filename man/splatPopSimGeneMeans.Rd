% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splatPop-simulate.R
\name{splatPopSimGeneMeans}
\alias{splatPopSimGeneMeans}
\title{Simulate gene means for splatPop}
\usage{
splatPopSimGeneMeans(sim, params, base.means.gene)
}
\arguments{
\item{sim}{SingleCellExperiment to add gene means to.}

\item{params}{SplatParams object with simulation parameters.}

\item{base.means.gene}{List of gene means for sample from matrix
generated by `splatPopSimulateMeans` and with the sample specified
in `splatPopSimulateSC`.}
}
\value{
SingleCellExperiment with simulated gene means.
}
\description{
Simulate outlier expression factors for splatPop. Genes with an outlier
factor not equal to 1 are replaced with the median mean expression
multiplied by the outlier factor.
}
\keyword{internal}
