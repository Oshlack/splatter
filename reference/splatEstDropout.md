# Estimate Splat dropout parameters

Estimate the midpoint and shape parameters for the logistic function
used when simulating dropout.

## Usage

``` r
splatEstDropout(norm.counts, params)
```

## Arguments

- norm.counts:

  library size normalised counts matrix.

- params:

  SplatParams object to store estimated values in.

## Value

SplatParams object with estimated values.

## Details

Logistic function parameters are estimated by fitting a logistic
function to the relationship between log2 mean gene expression and the
proportion of zeros in each gene. See
[`nls`](https://rdrr.io/r/stats/nls.html) for details of fitting. Note
this is done on the experiment level, more granular (eg. group or cell)
level dropout is not estimated.
