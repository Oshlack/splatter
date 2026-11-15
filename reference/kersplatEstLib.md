# Estimate Kersplat library size parameters

Estimate the library size parameters for the Kersplat simulation

## Usage

``` r
kersplatEstLib(counts, params, verbose)
```

## Arguments

- counts:

  counts matrix.

- params:

  KersplatParams object to store estimated values in.

- verbose:

  logical. Whether to print progress messages

## Value

KersplatParams object with library size parameters

## Details

Parameters for the log-normal distribution are estimated by fitting the
library sizes using
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).
All the fitting methods are tried and the fit with the best Cramer-von
Mises statistic is selected. The density of the library sizes is also
estimated using [`density`](https://rdrr.io/r/stats/density.html).
