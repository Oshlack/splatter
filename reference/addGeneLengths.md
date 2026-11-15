# Add gene lengths

Add gene lengths to an SingleCellExperiment object

## Usage

``` r
addGeneLengths(
  sce,
  method = c("generate", "sample"),
  loc = 7.9,
  scale = 0.7,
  lengths = NULL
)
```

## Arguments

- sce:

  SingleCellExperiment to add gene lengths to.

- method:

  Method to use for creating lengths.

- loc:

  Location parameter for the generate method.

- scale:

  Scale parameter for the generate method.

- lengths:

  Vector of lengths for the sample method.

## Value

SingleCellExperiment with added gene lengths

## Details

This function adds simulated gene lengths to the
[`rowData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
slot of a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object that can be used for calculating length normalised expression
values such as TPM or FPKM. The `generate` method simulates lengths
using a (rounded) log-normal distribution, with the default `loc` and
`scale` parameters based on human protein-coding genes. Alternatively
the `sample` method can be used which randomly samples lengths (with
replacement) from a supplied vector.

## Examples

``` r
# Default generate method
sce <- simpleSimulate()
#> Simulating means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.65 * dense matrix
sce <- addGeneLengths(sce)
head(rowData(sce))
#> DataFrame with 6 rows and 3 columns
#>              Gene    GeneMean    Length
#>       <character>   <numeric> <numeric>
#> Gene1       Gene1 7.17682e-05      1691
#> Gene2       Gene2 4.64094e+00      6187
#> Gene3       Gene3 1.35683e+00      3583
#> Gene4       Gene4 5.89116e-02      4134
#> Gene5       Gene5 5.62611e+00      3019
#> Gene6       Gene6 6.26506e-01      5804
# Sample method (human coding genes)
if (FALSE) { # \dontrun{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx.lens <- transcriptLengths(txdb, with.cds_len = TRUE)
tx.lens <- tx.lens[tx.lens$cds_len > 0, ]
gene.lens <- max(splitAsList(tx.lens$tx_len, tx.lens$gene_id))
sce <- addGeneLengths(sce, method = "sample", lengths = gene.lens)
} # }
```
