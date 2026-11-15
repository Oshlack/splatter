# Generate mock gff

Quick function to generate a mock gff.

## Usage

``` r
mockGFF(n.genes = 50, chromosome = 1, chr.length = 2e+06, seed = NULL)
```

## Arguments

- n.genes:

  Number of genes in mock gff file

- chromosome:

  Chromosome name

- chr.length:

  Length of mock chromosome

- seed:

  Optional: seed for random seed

## Value

data.frame containing mock gff data.

## Examples

``` r
gff <- mockGFF()
```
