% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-data.R
\name{mockBulkMatrix}
\alias{mockBulkMatrix}
\title{Generate mock bulk population scale expression data}
\usage{
mockBulkMatrix(n.genes = 100, n.samples = 50, seed = NULL)
}
\arguments{
\item{n.genes}{Number of genes in mock bulk data.}

\item{n.samples}{Number of samples in mock bulk data.}

\item{seed}{Optional: seed for random seed}
}
\value{
matrix containing mock bulk expression data.
}
\description{
Quick function to generate mock bulk expression data for a population, with
parameters estimated using real thyroid tissue data from GTEx.
}
\examples{
bulk <- mockBulkMatrix

}
