# Estimate gene mean and gene mean variance parameters

Estimate gene mean and gene mean variance parameters

## Usage

``` r
splatPopEstimateMeanCV(params, emp.gene.means)
```

## Arguments

- params:

  SplatPopParams object containing parameters for the simulation of the
  mean expression levels for the population. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- emp.gene.means:

  data.frame of empirical gene means across a population, where rows are
  genes and columns are individuals.

## Value

params object with estimated values.

## Details

Parameters for the mean gamma distribution are estimated by fitting the
mean (across the population) expression of genes that meet the criteria
(\<50% of samples have exp \<0.1) and parameters for the cv gamma
distribution are estimated for each bin of mean expression using the cv
of expression across the population for genes in that bin. Both are fit
using
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).
The "Nelder-Mead" method is used to fit the mean gamma distribution and
the maximum goodness-of-fit estimation method is used to minimise the
Cramer-von Mises distance for the CV distribution.
