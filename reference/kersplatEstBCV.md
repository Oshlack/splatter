# Estimate Kersplat BCV parameters

Estimate Biological Coefficient of Variation (BCV) parameters for the
Kersplat simulation

## Usage

``` r
kersplatEstBCV(counts, params, verbose)
```

## Arguments

- counts:

  counts matrix.

- params:

  KersplatParams object to store estimated values in.

- verbose:

  logical. Whether to print progress messages

## Value

KersplatParams object with estimated BCV parameters

## Details

The [`estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)
function is used to estimate the common dispersion across the dataset.
An exponential correction is applied based on fitting an exponential
relationship between simulated and estimated values. If this results in
a negative dispersion a simpler linear correction is applied instead.
