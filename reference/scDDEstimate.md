# Estimate scDD simulation parameters

Estimate simulation parameters for the scDD simulation from a real
dataset.

## Usage

``` r
scDDEstimate(
  counts,
  params = newSCDDParams(),
  verbose = TRUE,
  BPPARAM = SerialParam(),
  ...
)

# S3 method for class 'matrix'
scDDEstimate(
  counts,
  params = newSCDDParams(),
  verbose = TRUE,
  BPPARAM = SerialParam(),
  conditions,
  ...
)

# S3 method for class 'SingleCellExperiment'
scDDEstimate(
  counts,
  params = newSCDDParams(),
  verbose = TRUE,
  BPPARAM = SerialParam(),
  condition = "condition",
  ...
)

# Default S3 method
scDDEstimate(
  counts,
  params = newSCDDParams(),
  verbose = TRUE,
  BPPARAM = SerialParam(),
  condition,
  ...
)
```

## Arguments

- counts:

  either a counts matrix or a SingleCellExperiment object containing
  count data to estimate parameters from.

- params:

  SCDDParams object to store estimated values in.

- verbose:

  logical. Whether to show progress messages.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance giving the parallel back-end to be used. Default is
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  which uses a single core.

- ...:

  further arguments passed to or from other methods.

- conditions:

  Vector giving the condition that each cell belongs to. Conditions can
  be 1 or 2.

- condition:

  String giving the column that represents biological group of interest.

## Value

SCDDParams object containing the estimated parameters.

## Details

This function applies
[`preprocess`](https://rdrr.io/pkg/scDD/man/preprocess.html) to the
counts then uses [`scDD`](https://rdrr.io/pkg/scDD/man/scDD.html) to
estimate the numbers of each gene type to simulate. The output is then
converted to a SCDDParams object. See
[`preprocess`](https://rdrr.io/pkg/scDD/man/preprocess.html) and
[`scDD`](https://rdrr.io/pkg/scDD/man/scDD.html) for details.

## Examples

``` r
if (requireNamespace("scDD", quietly = TRUE)) {
    library(scuttle)
    set.seed(1)
    sce <- mockSCE(ncells = 20, ngenes = 100)

    colData(sce)$condition <- sample(1:2, ncol(sce), replace = TRUE)
    params <- scDDEstimate(sce, condition = "condition")
    params
}
#> Performing Median Normalization
#> Notice: 10 genes have less than 3 nonzero cells per condition.  Skipping these genes.
#> Setting up parallel back-end using 1 cores
#> Clustering observed expression data for each gene
#> Notice: Number of permutations is set to zero; using 
#>             Kolmogorov-Smirnov to test for differences in distributions
#>             instead of the Bayes Factor permutation test
#> Classifying significant genes into patterns
#> A Params object of class SCDDParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>     100       10   489515 
#> 
#> 11 additional parameters 
#> 
#> Data: 
#> (SCdat) 
#> SingleCellExperiment with 95 features and 20 cells 
#> 
#> Genes: 
#> (NDE)  (NDP)  (NDM)  (NDP)  (NEE)  (NEP) 
#>     0      0      0      0     95      5 
#> 
#> Fold change: 
#> [SD Range]   [Mode FC] 
#>       1, 3     2, 3, 4 
#> 
#> Variance: 
#> [Inflation] 
#>        1, 1 
#> 
#> Condition: 
#> [Condition] 
#>   condition 
#> 
```
