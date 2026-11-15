# Simulate group differential expression

Simulate differential expression. Differential expression factors for
each group are produced using
[`getLNormFactors`](http://oshlacklab.com/splatter/reference/getLNormFactors.md)
and these are added along with updated means for each group. For paths
care is taken to make sure they are simulated in the correct order.

## Usage

``` r
splatSimGroupDE(sim, params)

splatSimPathDE(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add differential expression to.

- params:

  splatParams object with simulation parameters.

## Value

SingleCellExperiment with simulated differential expression.
