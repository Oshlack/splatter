# Generate mock eQTL mapping results

Quick function to generate mock eQTL mapping results, with parameters
estimated using real eQTL mapping results from GTEx using thyroid
tissue.

## Usage

``` r
mockBulkeQTL(n.genes = 500, seed = NULL)
```

## Arguments

- n.genes:

  Number of genes in mock eQTL data.

- seed:

  Optional: seed for random seed

## Value

data.frame containing mock bulk eQTL mapping results.

## Examples

``` r
eqtl <- mockBulkeQTL()
```
