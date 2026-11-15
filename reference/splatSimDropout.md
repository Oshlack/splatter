# Simulate dropout

A logistic function is used to form a relationship between the
expression level of a gene and the probability of dropout, giving a
probability for each gene in each cell. These probabilities are used in
a Bernoulli distribution to decide which counts should be dropped.

## Usage

``` r
splatSimDropout(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add dropout to.

- params:

  SplatParams object with simulation parameters.

## Value

SingleCellExperiment with simulated dropout and observed counts.
