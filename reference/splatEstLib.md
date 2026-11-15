# Estimate Splat library size parameters

The Shapiro-Wilks test is used to determine if the library sizes are
normally distributed. If so a normal distribution is fitted to the
library sizes, if not (most cases) a log-normal distribution is fitted
and the estimated parameters are added to the params object. See
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
for details on the fitting.

## Usage

``` r
splatEstLib(counts, params)
```

## Arguments

- counts:

  counts matrix to estimate parameters from.

- params:

  splatParams object to store estimated values in.

## Value

SplatParams object with estimated values.
