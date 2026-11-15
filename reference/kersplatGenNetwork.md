# Generate Kersplat gene network

Generate a gene network for the Kersplat simulation

## Usage

``` r
kersplatGenNetwork(params, verbose)
```

## Arguments

- params:

  KersplatParams object containing simulation parameters.

- verbose:

  logical. Whether to print progress messages

## Value

KersplatParams object with gene network

## Details

Currently a very simple approach is used which needs to be improved. A
network is generated using the
[`sample_forestfire`](https://r.igraph.org/reference/sample_forestfire.html)
function and edge weights are sampled from a standard normal
distribution.
