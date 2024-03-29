% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kersplat-estimate.R
\name{kersplatEstLib}
\alias{kersplatEstLib}
\title{Estimate Kersplat library size parameters}
\usage{
kersplatEstLib(counts, params, verbose)
}
\arguments{
\item{counts}{counts matrix.}

\item{params}{KersplatParams object to store estimated values in.}

\item{verbose}{logical. Whether to print progress messages}
}
\value{
KersplatParams object with library size parameters
}
\description{
Estimate the library size parameters for the Kersplat simulation
}
\details{
Parameters for the log-normal distribution are estimated by fitting the
library sizes using \code{\link[fitdistrplus]{fitdist}}. All the fitting
methods are tried and the fit with the best Cramer-von Mises statistic is
selected. The density of the library sizes is also estimated using
\code{\link[stats]{density}}.
}
\keyword{internal}
