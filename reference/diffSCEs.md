# Diff SingleCellExperiment objects

Combine the data from several SingleCellExperiment objects and produce
some basic plots comparing them to a reference.

## Usage

``` r
diffSCEs(
  sces,
  ref,
  point.size = 0.1,
  point.alpha = 0.1,
  fits = TRUE,
  colours = NULL
)
```

## Arguments

- sces:

  named list of SingleCellExperiment objects to combine and compare.

- ref:

  string giving the name of the SingleCellExperiment to use as the
  reference

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

This function aims to look at the differences between a reference
SingleCellExperiment and one or more others. It requires each
SingleCellExperiment to have the same dimensions. Properties are
compared by ranks, for example when comparing the means the values are
ordered and the differences between the reference and another dataset
plotted. A series of Q-Q plots are also returned.

The returned list has five items:

- `Reference`:

  The SingleCellExperiment used as the reference.

- `RowData`:

  Combined feature data from the provided SingleCellExperiments.

- `ColData`:

  Combined column data from the provided SingleCellExperiments.

- `Plots`:

  Difference plots

  `Means`

  : Boxplot of mean differences.

  `Variances`

  : Boxplot of variance differences.

  `MeanVar`

  : Scatter plot showing the difference from the reference variance
    across expression ranks.

  `LibraeySizes`

  : Boxplot of the library size differences.

  `ZerosGene`

  : Boxplot of the differences in the percentage of each gene that is
    zero.

  `ZerosCell`

  : Boxplot of the differences in the percentage of each cell that is
    zero.

  `MeanZeros`

  : Scatter plot showing the difference from the reference percentage of
    zeros across expression ranks.

- `QQPlots`:

  Quantile-Quantile plots

  `Means`

  : Q-Q plot of the means.

  `Variances`

  : Q-Q plot of the variances.

  `LibrarySizes`

  : Q-Q plot of the library sizes.

  `ZerosGene`

  : Q-Q plot of the percentage of zeros per gene.

  `ZerosCell`

  : Q-Q plot of the percentage of zeros per cell.

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
#> Skipping 'TrueCounts': estimated sparse size 2.6 * dense matrix
#> Skipping 'counts': estimated sparse size 2.6 * dense matrix
#> Done!
sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#> Simulating means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.63 * dense matrix
difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")
names(difference)
#> [1] "Reference" "RowData"   "ColData"   "Plots"     "QQPlots"  
names(difference$Plots)
#> [1] "Means"        "Variances"    "MeanVar"      "LibrarySizes" "ZerosGene"   
#> [6] "ZerosCell"    "MeanZeros"   
```
