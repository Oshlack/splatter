% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splat-estimate.R
\name{splatEstBCV}
\alias{splatEstBCV}
\title{Estimate Splat Biological Coefficient of Variation parameters}
\usage{
splatEstBCV(counts, params)
}
\arguments{
\item{counts}{counts matrix to estimate parameters from.}

\item{params}{SplatParams object to store estimated values in.}
}
\value{
SplatParams object with estimated values.
}
\description{
Parameters are estimated using the \code{\link[edgeR]{estimateDisp}} function
in the \code{edgeR} package.
}
\details{
The \code{\link[edgeR]{estimateDisp}} function is used to estimate the common
dispersion and prior degrees of freedom. See
\code{\link[edgeR]{estimateDisp}} for details. When estimating parameters on
simulated data we found a broadly linear relationship between the true
underlying common dispersion and the \code{edgR} estimate, therefore we
apply a small correction, \code{disp = 0.1 + 0.25 * edgeR.disp}.
}
\keyword{internal}
