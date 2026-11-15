# Make difference panel

Combine the plots from `diffSCEs` into a single panel.

## Usage

``` r
makeDiffPanel(
  diff,
  title = "Difference comparison",
  labels = c("Means", "Variance", "Library size", "Zeros per cell", "Zeros per gene",
    "Mean-variance relationship", "Mean-zeros relationship")
)
```

## Arguments

- diff:

  list returned by
  [`diffSCEs`](http://oshlacklab.com/splatter/reference/diffSCEs.md).

- title:

  title for the panel.

- labels:

  vector of labels for each of the seven sections.

## Value

Combined panel plot

## Examples

``` r
# \donttest{
sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#> Getting parameters...
#> Creating simulation object...
#> Simulating library sizes...
#> Simulating gene means...
#> Simulating BCV...
#> Simulating counts...
#> Simulating dropout (if needed)...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BCV': estimated sparse size 1.5 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.49 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 2.6 * dense matrix
#> Skipping 'counts': estimated sparse size 2.6 * dense matrix
#> Done!
sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#> Simulating means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.65 * dense matrix
difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")
panel <- makeDiffPanel(difference)
#> Notch went outside hinges
#> ℹ Do you want `notch = FALSE`?
# }
```
