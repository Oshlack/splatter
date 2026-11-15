# Estimate population/eQTL simulation parameters

Estimate simulation parameters for the eQTL population simulation from
real data. See the individual estimation functions for more details on
how this is done.

## Usage

``` r
splatPopEstimate(
  counts = NULL,
  means = NULL,
  eqtl = NULL,
  params = newSplatPopParams()
)
```

## Arguments

- counts:

  either a counts matrix or a SingleCellExperiment object containing
  count data to estimate parameters from.

- means:

  Matrix of real gene means across a population, where each row is a
  gene and each column is an individual in the population.

- eqtl:

  data.frame with all or top eQTL pairs from a real eQTL analysis. Must
  include columns: 'gene_id', 'pval_nominal', and 'slope'.

- params:

  SplatPopParams object containing parameters for the simulation of the
  mean expression levels for the population. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

## Value

SplatPopParams object containing the estimated parameters.

## See also

[`splatPopEstimateEffectSize`](http://oshlacklab.com/splatter/reference/splatPopEstimateEffectSize.md),
[`splatPopEstimateMeanCV`](http://oshlacklab.com/splatter/reference/splatPopEstimateMeanCV.md)

## Examples

``` r

if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {
    # Load example data
    library(scuttle)

    sce <- mockSCE()
    params <- splatPopEstimate(sce)
}
```
