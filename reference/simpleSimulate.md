# Simple simulation

Simulate counts from a simple negative binomial distribution without
simulated library sizes, differential expression etc.

## Usage

``` r
simpleSimulate(
  params = newSimpleParams(),
  sparsify = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- params:

  SimpleParams object containing simulation parameters.

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment containing simulated counts

## Details

Gene means are simulated from a gamma distribution with
`shape = mean.shape` and `rate = mean.rate`. Counts are then simulated
from a negative binomial distribution with `mu = means` and
`size = 1 / counts.disp`. See
[`SimpleParams`](http://oshlacklab.com/splatter/reference/SimpleParams.md)
for more details of the parameters.

## Examples

``` r
sim <- simpleSimulate()
#> Simulating means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.65 * dense matrix
# Override default parameters
sim <- simpleSimulate(nGenes = 1000, nCells = 50)
#> Simulating means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.66 * dense matrix
```
