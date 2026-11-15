# splatPopParseEmpirical

Parse splatPop key information from empirical data provided.

## Usage

``` r
splatPopParseEmpirical(
  vcf = vcf,
  gff = gff,
  eqtl = eqtl,
  means = means,
  params = params
)
```

## Arguments

- vcf:

  VariantAnnotation object containing genotypes of samples.

- gff:

  Either NULL or a data.frame object containing a GFF/GTF file.

- eqtl:

  Either NULL or if simulating population parameters directly from
  empirical data, a data.frame with empirical/desired eQTL results. To
  see required format, run
  [`mockEmpiricalSet()`](http://oshlacklab.com/splatter/reference/mockEmpiricalSet.md)
  and see eqtl output.

- means:

  Either NULL or if simulating population parameters directly from
  empirical data, a Matrix of real gene means across a population, where
  each row is a gene and each column is an individual in the population.
  To see required format, run
  [`mockEmpiricalSet()`](http://oshlacklab.com/splatter/reference/mockEmpiricalSet.md)
  and see means output.

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

## Value

A partial splatPop `key`

## Details

NOTE: This function will cause some of the parameters in the
splatPopParams object to be ignored, such as population level gene mean
and variance and eQTL parameters.

This function will ignore a number of parameters defined in
splatPopParams, instead pulling key information directly from provided
VCF, GFF, gene means, and eQTL mapping result data provided.
