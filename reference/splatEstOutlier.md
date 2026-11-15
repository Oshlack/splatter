# Estimate Splat expression outlier parameters

Parameters are estimated by comparing means of individual genes to the
median mean expression level.

## Usage

``` r
splatEstOutlier(norm.counts, params)
```

## Arguments

- norm.counts:

  library size normalised counts matrix.

- params:

  SplatParams object to store estimated values in.

## Value

SplatParams object with estimated values.

## Details

Expression outlier genes are detected using the Median Absolute
Deviation (MAD) from median method. If the log2 mean expression of a
gene is greater than two MADs above the median log2 mean expression it
is designated as an outlier. The proportion of outlier genes is used to
estimate the outlier probability. Factors for each outlier gene are
calculated by dividing mean expression by the median mean expression. A
log-normal distribution is then fitted to these factors in order to
estimate the outlier factor location and scale parameters using
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).
