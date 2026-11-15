# Simulate batch means

Simulate a mean for each gene in each cell incorporating batch effect
factors.

## Usage

``` r
splatSimBatchCellMeans(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add batch means to.

- params:

  SplatParams object with simulation parameters.

## Value

SingleCellExperiment with simulated batch means.
