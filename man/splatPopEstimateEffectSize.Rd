% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splatPop-estimate.R
\name{splatPopEstimateEffectSize}
\alias{splatPopEstimateEffectSize}
\title{Estimate eQTL Effect Size parameters}
\usage{
splatPopEstimateEffectSize(params, eqtl)
}
\arguments{
\item{params}{SplatPopParams object containing parameters for the
simulation of the mean expression levels for the population.
See \code{\link{SplatPopParams}} for details.}

\item{eqtl}{data.frame with all or top eQTL pairs from a real eQTL analysis.
Must include columns: gene_id, pval_nominal, and slope.}
}
\value{
params object with estimated values.
}
\description{
Estimate rate and shape parameters for the gamma distribution used to
simulate eQTL (eSNP-eGene) effect sizes.
}
\details{
Parameters for the gamma distribution are estimated by fitting the top eSNP-
eGene pair effect sizes using \code{\link[fitdistrplus]{fitdist}}. The
maximum goodness-of-fit estimation method is used to minimise the
Cramer-von Mises distance. This can fail in some situations, in which case
the method of moments estimation method is used instead.
}
\keyword{internal}
