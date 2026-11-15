# Compare SingleCellExperiment objects

Combine the data from several SingleCellExperiment objects and produce
some basic plots comparing them.

## Usage

``` r
compareSCEs(
  sces,
  point.size = 0.1,
  point.alpha = 0.1,
  fits = TRUE,
  colours = NULL
)
```

## Arguments

- sces:

  named list of SingleCellExperiment objects to combine and compare.

- point.size:

  size of points in scatter plots.

- point.alpha:

  opacity of points in scatter plots.

- fits:

  whether to include fits in scatter plots.

- colours:

  vector of colours to use for each dataset.

## Value

List containing the combined datasets and plots.

## Details

The returned list has three items:

- `RowData`:

  Combined row data from the provided SingleCellExperiments.

- `ColData`:

  Combined column data from the provided SingleCellExperiments.

- `Plots`:

  Comparison plots

  `Means`

  : Boxplot of mean distribution.

  `Variances`

  : Boxplot of variance distribution.

  `MeanVar`

  : Scatter plot with fitted lines showing the mean-variance
    relationship.

  `LibrarySizes`

  : Boxplot of the library size distribution.

  `ZerosGene`

  : Boxplot of the percentage of each gene that is zero.

  `ZerosCell`

  : Boxplot of the percentage of each cell that is zero.

  `MeanZeros`

  : Scatter plot with fitted lines showing the mean-zeros relationship.

  `VarGeneCor`

  : Heatmap of correlation of the 100 most variable genes.

The plots returned by this function are created using
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) and are
only a sample of the kind of plots you might like to consider. The data
used to create these plots is also returned and should be in the correct
format to allow you to create further plots using
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Examples

``` r
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
#> Skipping 'TrueCounts': estimated sparse size 2.59 * dense matrix
#> Skipping 'counts': estimated sparse size 2.59 * dense matrix
#> Done!
sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#> Simulating means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.6 * dense matrix
comparison <- compareSCEs(list(Splat = sim1, Simple = sim2))
names(comparison)
#> [1] "RowData" "ColData" "Plots"  
names(comparison$Plots)
#> [1] "Means"        "Variances"    "MeanVar"      "LibrarySizes" "ZerosGene"   
#> [6] "ZerosCell"    "MeanZeros"    "VarGeneCor"  
```
