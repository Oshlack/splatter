# Simulate BCV means

Simulate means for each gene in each cell that are adjusted to follow a
mean-variance trend using Biological Coefficient of Variation taken from
and inverse gamma distribution.

## Usage

``` r
splatSimBCVMeans(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add BCV means to.

- params:

  SplatParams object with simulation parameters.

## Value

SingleCellExperiment with simulated BCV means.
