# Simulate mean gene expression matrix without eQTL effects

Gene mean expression levels are assigned to each gene for each pair
randomly from a normal distribution parameterized using the mean and cv
assigned to each gene in the key. If gene means matrix is provided,
those will be used instead.

## Usage

``` r
splatPopSimMeans(vcf, key, means)
```

## Arguments

- vcf:

  VariantAnnotation object containing genotypes of samples.

- key:

  Partial splatPop key data.frame.

- means:

  Null or matrix of gene means to use

## Value

matrix of gene mean expression levels WITHOUT eQTL effects.
