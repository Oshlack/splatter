# splatPopSimulateMeans

Simulate mean expression levels for all genes for all samples, with
between sample correlation structure simulated with eQTL effects and
with the option to simulate multiple groups (i.e. cell-types).

## Usage

``` r
splatPopSimulateMeans(
  vcf = mockVCF(),
  params = newSplatPopParams(nGenes = 1000),
  verbose = TRUE,
  key = NULL,
  gff = NULL,
  eqtl = NULL,
  means = NULL,
  ...
)
```

## Arguments

- vcf:

  VariantAnnotation object containing genotypes of samples.

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- verbose:

  logical. Whether to print progress messages.

- key:

  Either FALSE or a data.frame object containing a full or partial
  splatPop key.

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

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

A list containing: `means` a matrix (or list of matrices if n.groups
\> 1) with the simulated mean gene expression value for each gene (row)
and each sample (column), `key` a data.frame with population information
including eQTL and group effects, and `condition` a named array
containing conditional group assignments for each sample.

## Details

SplatPopParams can be set in a variety of ways. 1. If not provided,
default parameters are used. 2. Default parameters can be overridden by
supplying desired parameters using
[`setParams`](http://oshlacklab.com/splatter/reference/setParams.md). 3.
Parameters can be estimated from real data of your choice using
[`splatPopEstimate`](http://oshlacklab.com/splatter/reference/splatPopEstimate.md).

`splatPopSimulateMeans` involves the following steps:

1.  Load population key or generate random or GFF/GTF based key.

2.  Format and subset genotype data from the VCF file.

3.  If not in key, assign expression mean and variance to each gene.

4.  If not in key, assign eGenes-eSNPs pairs and effect sizes.

5.  If not in key and groups \>1, assign subset of eQTL associations as
    group-specific and assign DEG group effects.

6.  Simulate mean gene expression matrix without eQTL effects

7.  Quantile normalize by sample to fit single-cell expression
    distribution as defined in `splatEstimate`.

8.  Add quantile normalized gene mean and cv info the eQTL key.

9.  Add eQTL effects to means matrix.

## See also

[`splatPopParseVCF`](http://oshlacklab.com/splatter/reference/splatPopParseVCF.md),
[`splatPopParseGenes`](http://oshlacklab.com/splatter/reference/splatPopParseGenes.md),
[`splatPopAssignMeans`](http://oshlacklab.com/splatter/reference/splatPopAssignMeans.md),
[`splatPopQuantNorm`](http://oshlacklab.com/splatter/reference/splatPopQuantNorm.md),
[`splatPopQuantNormKey`](http://oshlacklab.com/splatter/reference/splatPopQuantNormKey.md)
[`splatPopeQTLEffects`](http://oshlacklab.com/splatter/reference/splatPopeQTLEffects.md),
[`splatPopGroupEffects`](http://oshlacklab.com/splatter/reference/splatPopGroupEffects.md),
[`splatPopSimMeans`](http://oshlacklab.com/splatter/reference/splatPopSimMeans.md),
[`splatPopSimEffects`](http://oshlacklab.com/splatter/reference/splatPopSimEffects.md),

## Examples

``` r
# \donttest{
if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {
    means <- splatPopSimulateMeans()
}
#> Simulating data for genes in GFF...
#> Simulating gene means for population...
# }
```
