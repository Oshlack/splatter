# Generate mock vcf

Quick function to generate mock vcf file. Note this data has unrealistic
population structure.

## Usage

``` r
mockVCF(
  n.snps = 200,
  n.samples = 5,
  chromosome = 1,
  chr.length = 2e+06,
  seed = NULL
)
```

## Arguments

- n.snps:

  Number of SNPs in mock vcf file.

- n.samples:

  Number of samples in mock bulk data.

- chromosome:

  Chromosome name

- chr.length:

  Length of mock chromosome

- seed:

  Optional: seed for random seed

## Value

data.frame containing mock vcf data.

## Examples

``` r
vcf <- mockVCF()
```
