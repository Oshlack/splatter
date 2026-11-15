# Format and subset genotype data from a VCF file.

Extract numeric alleles from vcf object and filter out SNPs missing
genotype data or outside the Minor Allele Frequency range in
`SplatPopParams`.

## Usage

``` r
splatPopParseVCF(vcf, params)
```

## Arguments

- vcf:

  VariantAnnotation object containing genotypes of samples.

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

## Value

Genotype data.frame
