# Estimate PhenoPath simulation parameters

Estimate simulation parameters for the PhenoPath simulation from a real
dataset.

## Usage

``` r
phenoEstimate(counts, params = newPhenoParams())

# S3 method for class 'SingleCellExperiment'
phenoEstimate(counts, params = newPhenoParams())

# S3 method for class 'matrix'
phenoEstimate(counts, params = newPhenoParams())
```

## Arguments

- counts:

  either a counts matrix or an SingleCellExperiment object containing
  count data to estimate parameters from.

- params:

  PhenoParams object to store estimated values in.

## Value

PhenoParams object containing the estimated parameters.

## Details

The `nGenes` and `nCells` parameters are taken from the size of the
input data. The total number of genes is evenly divided into the four
types. See
[`PhenoParams`](http://oshlacklab.com/splatter/reference/PhenoParams.md)
for more details on the parameters.

## Examples

``` r
if (requireNamespace("phenopath", quietly = TRUE)) {
    # Load example data
    library(scuttle)
    set.seed(1)
    sce <- mockSCE()

    params <- phenoEstimate(sce)
    params
}
#> A Params object of class PhenoParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>    2000      200   787110 
#> 
#> 4 additional parameters 
#> 
#> Genes: 
#>              [DE]              [PST]       [PST + BETA]  [DE + PST + BETA] 
#>               500                500                500                500 
#> 
```
