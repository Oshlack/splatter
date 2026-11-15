# Simulate batch effects

Simulate batch effects. Batch effect factors for each batch are produced
using
[`getLNormFactors`](http://oshlacklab.com/splatter/reference/getLNormFactors.md)
and these are added along with updated means for each batch.

## Usage

``` r
splatPopSimBatchEffects(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add batch effects to.

- params:

  SplatParams object with simulation parameters.

## Value

SingleCellExperiment with simulated batch effects.
