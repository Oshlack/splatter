# Summarise diffSCEs

Summarise the results of
[`diffSCEs`](http://oshlacklab.com/splatter/reference/diffSCEs.md).
Calculates the Median Absolute Deviation (MAD), Mean Absolute Error
(MAE), Root Mean Squared Error (RMSE) and Kolmogorov-Smirnov (KS)
statistics for the various properties and ranks them.

## Usage

``` r
summariseDiff(diff)
```

## Arguments

- diff:

  Output from
  [`diffSCEs`](http://oshlacklab.com/splatter/reference/diffSCEs.md)

## Value

data.frame with MADs, MAEs, RMSEs, scaled statistics and ranks

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
#> Converting 'counts' to sparse matrix: estimated sparse size 0.66 * dense matrix
difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")
summary <- summariseDiff(difference)
head(summary)
#>   Dataset Statistic          MAD MADScaled MADRank         MAE MAEScaled
#> 1   Splat      Mean     2.364391       NaN       1     2.46779       NaN
#> 2   Splat  Variance    11.503897       NaN       1    10.28866       NaN
#> 3   Splat ZerosGene    35.000000       NaN       1    41.65500       NaN
#> 4   Splat   MeanVar    10.759043       NaN       1    11.85966       NaN
#> 5   Splat MeanZeros    42.500000       NaN       1    42.40500       NaN
#> 6   Splat   LibSize 61409.000000       NaN       1 61161.40000       NaN
#>   MAERank         RMSE RMSEScaled RMSERank    KS        KSPVal KSRank
#> 1       1     2.951595        NaN        1 0.339  2.463053e-50      1
#> 2       1    12.857164        NaN        1 0.592 1.248673e-152      1
#> 3       1    45.557381        NaN        1 0.561 4.163287e-137      1
#> 4       1    14.895204        NaN        1    NA            NA     NA
#> 5       1    51.702756        NaN        1    NA            NA     NA
#> 6       1 62260.204076        NaN        1 1.000  1.450889e-11      1
```
