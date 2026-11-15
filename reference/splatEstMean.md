# Estimate Splat mean parameters

Estimate rate and shape parameters for the gamma distribution used to
simulate gene expression means.

## Usage

``` r
splatEstMean(norm.counts, params)
```

## Arguments

- norm.counts:

  library size normalised counts matrix.

- params:

  SplatParams object to store estimated values in.

## Value

SplatParams object containing the estimated parameters.

## Details

Parameters for the gamma distribution are estimated by fitting the mean
normalised counts using
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).
The 'maximum goodness-of-fit estimation' method is used to minimise the
Cramer-von Mises distance. This can fail in some situations, in which
case the 'method of moments estimation' method is used instead. Prior to
fitting the means are winsorized by setting the top and bottom 10
percent of values to the 10th and 90th percentiles.
