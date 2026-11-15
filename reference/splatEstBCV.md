# Estimate Splat Biological Coefficient of Variation parameters

Parameters are estimated using the
[`estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)
function in the `edgeR` package.

## Usage

``` r
splatEstBCV(counts, params)
```

## Arguments

- counts:

  counts matrix to estimate parameters from.

- params:

  SplatParams object to store estimated values in.

## Value

SplatParams object with estimated values.

## Details

The [`estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)
function is used to estimate the common dispersion and prior degrees of
freedom. See
[`estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html) for
details. When estimating parameters on simulated data we found a broadly
linear relationship between the true underlying common dispersion and
the `edgR` estimate, therefore we apply a small correction,
`disp = 0.1 + 0.25 * edgeR.disp`.
