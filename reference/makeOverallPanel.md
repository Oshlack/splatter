# Make overall panel

Combine the plots from `compSCEs` and `diffSCEs` into a single panel.

## Usage

``` r
makeOverallPanel(
  comp,
  diff,
  title = "Overall comparison",
  row.labels = c("Means", "Variance", "Mean-variance relationship", "Library size",
    "Zeros per cell", "Zeros per gene", "Mean-zeros relationship")
)
```

## Arguments

- comp:

  list returned by
  [`compareSCEs`](http://oshlacklab.com/splatter/reference/compareSCEs.md).

- diff:

  list returned by
  [`diffSCEs`](http://oshlacklab.com/splatter/reference/diffSCEs.md).

- title:

  title for the panel.

- row.labels:

  vector of labels for each of the seven rows.

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
#> Skipping 'TrueCounts': estimated sparse size 2.59 * dense matrix
#> Skipping 'counts': estimated sparse size 2.59 * dense matrix
#> Done!
sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#> Simulating means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.63 * dense matrix
comparison <- compareSCEs(list(Splat = sim1, Simple = sim2))
difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")
panel <- makeOverallPanel(comparison, difference)
#> Notch went outside hinges
#> ℹ Do you want `notch = FALSE`?
#> Notch went outside hinges
#> ℹ Do you want `notch = FALSE`?
#> Notch went outside hinges
#> ℹ Do you want `notch = FALSE`?
#> Warning: log-10 transformation introduced infinite values.
#> Warning: log-10 transformation introduced infinite values.
#> Warning: Removed 227 rows containing non-finite outside the scale range
#> (`stat_smooth()`).
# }
```
