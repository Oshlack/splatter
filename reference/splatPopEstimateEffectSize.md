# Estimate eQTL Effect Size parameters

Estimate rate and shape parameters for the gamma distribution used to
simulate eQTL (eSNP-eGene) effect sizes.

## Usage

``` r
splatPopEstimateEffectSize(params, eqtl)
```

## Arguments

- params:

  SplatPopParams object containing parameters for the simulation of the
  mean expression levels for the population. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- eqtl:

  data.frame with all or top eQTL pairs from a real eQTL analysis. Must
  include columns: gene_id, pval_nominal, and slope.

## Value

params object with estimated values.

## Details

Parameters for the gamma distribution are estimated by fitting the top
eSNP- eGene pair effect sizes using
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).
The maximum goodness-of-fit estimation method is used to minimise the
Cramer-von Mises distance. This can fail in some situations, in which
case the method of moments estimation method is used instead.
