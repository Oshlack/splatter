# Simulate gene means

Simulate gene means from a gamma distribution. Also simulates outlier
expression factors. Genes with an outlier factor not equal to 1 are
replaced with the median mean expression multiplied by the outlier
factor.

## Usage

``` r
splatSimGeneMeans(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add gene means to.

- params:

  SplatParams object with simulation parameters.

## Value

SingleCellExperiment with simulated gene means.
