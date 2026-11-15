# Generate mock bulk population scale expression data

Quick function to generate mock bulk expression data for a population,
with parameters estimated using real thyroid tissue data from GTEx.

## Usage

``` r
mockBulkMatrix(n.genes = 100, n.samples = 50, seed = NULL)
```

## Arguments

- n.genes:

  Number of genes in mock bulk data.

- n.samples:

  Number of samples in mock bulk data.

- seed:

  Optional: seed for random seed

## Value

matrix containing mock bulk expression data.

## Examples

``` r
bulk <- mockBulkMatrix
```
