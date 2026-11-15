# Quantile normalize by sample to fit sc expression distribution.

For each sample, expression values are quantile normalized (qgamma)
using the gamma distribution parameterized from splatEstimate(). This
ensures the simulated gene means reflect the distribution expected from
a sc dataset and not a bulk dataset.

## Usage

``` r
splatPopQuantNorm(params, means)
```

## Arguments

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- means:

  Mean gene expression matrix with eQTL effects.

## Value

matrix of quantile normalized gene mean expression levels.

## Examples

``` r

if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {
    bulk.means <- mockBulkMatrix(n.genes = 100, n.samples = 100)
    bulk.qnorm <- splatPopQuantNorm(newSplatPopParams(), bulk.means)
}
```
