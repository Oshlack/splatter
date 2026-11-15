# Simulate cell means

Simulate a gene by cell matrix giving the mean expression for each gene
in each cell. Cells start with the mean expression for the group they
belong to (when simulating groups) or cells are assigned the mean
expression from a random position on the appropriate path (when
simulating paths). The selected means are adjusted for each cell's
expected library size.

## Usage

``` r
splatSimSingleCellMeans(sim, params)

splatSimGroupCellMeans(sim, params)

splatSimPathCellMeans(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add cell means to.

- params:

  SplatParams object with simulation parameters.

## Value

SingleCellExperiment with added cell means.
