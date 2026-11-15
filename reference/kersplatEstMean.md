# Estimate Kersplat means

Estimate mean parameters for the Kersplat simulation

## Usage

``` r
kersplatEstMean(norm.counts, params, verbose)
```

## Arguments

- norm.counts:

  library size normalised counts matrix.

- params:

  KersplatParams object to store estimated values in.

- verbose:

  logical. Whether to print progress messages

## Value

KersplatParams object with estimated means

## Details

Parameters for the gamma distribution are estimated by fitting the mean
normalised counts using
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).
All the fitting methods are tried and the fit with the best Cramer-von
Mises statistic is selected. The density of the means is also estimated
using [`density`](https://rdrr.io/r/stats/density.html).

Expression outlier genes are detected using the Median Absolute
Deviation (MAD) from median method. If the log2 mean expression of a
gene is greater than two MADs above the median log2 mean expression it
is designated as an outlier. The proportion of outlier genes is used to
estimate the outlier probability. Factors for each outlier gene are
calculated by dividing mean expression by the median mean expression. A
log-normal distribution is then fitted to these factors in order to
estimate the outlier factor location and scale parameters using the
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
MLE method.
