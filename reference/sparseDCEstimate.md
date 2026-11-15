# Estimate SparseDC simulation parameters

Estimate simulation parameters for the SparseDC simulation from a real
dataset.

## Usage

``` r
sparseDCEstimate(
  counts,
  conditions,
  nclusters,
  norm = TRUE,
  params = newSparseDCParams()
)

# S3 method for class 'SingleCellExperiment'
sparseDCEstimate(
  counts,
  conditions,
  nclusters,
  norm = TRUE,
  params = newSparseDCParams()
)

# S3 method for class 'matrix'
sparseDCEstimate(
  counts,
  conditions,
  nclusters,
  norm = TRUE,
  params = newSparseDCParams()
)
```

## Arguments

- counts:

  either a counts matrix or an SingleCellExperiment object containing
  count data to estimate parameters from.

- conditions:

  numeric vector giving the condition each cell belongs to.

- nclusters:

  number of cluster present in the dataset.

- norm:

  logical, whether to library size normalise counts before estimation.
  Set this to FALSE if counts is already normalised.

- params:

  PhenoParams object to store estimated values in.

## Value

SparseParams object containing the estimated parameters.

## Details

The `nGenes` and `nCells` parameters are taken from the size of the
input data. The counts are preprocessed using
[`pre_proc_data`](https://rdrr.io/pkg/SparseDC/man/pre_proc_data.html)
and then parameters are estimated using
[`sparsedc_cluster`](https://rdrr.io/pkg/SparseDC/man/sparsedc_cluster.html)
using lambda values calculated using
[`lambda1_calculator`](https://rdrr.io/pkg/SparseDC/man/lambda1_calculator.html)
and
[`lambda2_calculator`](https://rdrr.io/pkg/SparseDC/man/lambda2_calculator.html).

See
[`SparseDCParams`](http://oshlacklab.com/splatter/reference/SparseDCParams.md)
for more details on the parameters.

## Examples

``` r
if (requireNamespace("SparseDC", quietly = TRUE)) {
    # Load example data
    library(scuttle)
    set.seed(1)
    sce <- mockSCE(ncells = 20, ngenes = 100)

    conditions <- sample(1:2, ncol(sce), replace = TRUE)

    params <- sparseDCEstimate(sce, conditions, nclusters = 3)
    params
}
#> A Params object of class SparseDCParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>     100       10   489515 
#> 
#> 7 additional parameters 
#> 
#> Markers: 
#> (NUMBER)  (SHARED)    [Same] 
#>        1         1     FALSE 
#> 
#> Clusters: 
#> (CONDITION 1)  (CONDITION 2) 
#>          2, 3        1, 2, 3 
#> 
#> Means: 
#> [Lower]  [Upper] 
#>       1        2 
#> 
```
