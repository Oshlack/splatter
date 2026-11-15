# Generate set of "empirical" mock data

Quick function to generate matching mock VCF, bulk expression, and eQTL
data, useful for running splatPopEmpiricalMeans

## Usage

``` r
mockEmpiricalSet(
  n.genes = 20,
  n.snps = 1000,
  n.samples = 10,
  chromosome = 1,
  chr.length = 2e+06,
  seed = NULL
)
```

## Arguments

- n.genes:

  Number of genes in mock eQTL data.

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

list(gff=mockGFF, vcf=mockVCF, means=mockMEANS, eqtl=mockEQTL)

## Examples

``` r
empirical <- mockEmpiricalSet()
```
