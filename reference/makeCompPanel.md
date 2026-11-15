# Make comparison panel

Combine the plots from `compareSCEs` into a single panel.

## Usage

``` r
makeCompPanel(
  comp,
  title = "Comparison",
  labels = c("Means", "Variance", "Mean-variance relationship", "Library size",
    "Zeros per gene", "Zeros per cell", "Mean-zeros relationship")
)
```

## Arguments

- comp:

  list returned by
  [`compareSCEs`](http://oshlacklab.com/splatter/reference/compareSCEs.md).

- title:

  title for the panel.

- labels:

  vector of labels for each of the seven plots.

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
#> Skipping 'CellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 2.66 * dense matrix
#> Skipping 'counts': estimated sparse size 2.66 * dense matrix
#> Done!
sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#> Simulating means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.63 * dense matrix
comparison <- compareSCEs(list(Splat = sim1, Simple = sim2))
panel <- makeCompPanel(comparison)
#> Notch went outside hinges
#> ℹ Do you want `notch = FALSE`?
#> Notch went outside hinges
#> ℹ Do you want `notch = FALSE`?
#> Notch went outside hinges
#> ℹ Do you want `notch = FALSE`?
#> Notch went outside hinges
#> ℹ Do you want `notch = FALSE`?
#> Warning: log-10 transformation introduced infinite values.
#> Warning: log-10 transformation introduced infinite values.
#> Warning: Removed 204 rows containing non-finite outside the scale range
#> (`stat_smooth()`).
# }
```
