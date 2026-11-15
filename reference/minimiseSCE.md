# Minimise SCE

Reduce the size of a SingleCellExperiment object by unneeded
information.

## Usage

``` r
minimiseSCE(
  sce,
  rowData.keep = FALSE,
  colData.keep = FALSE,
  metadata.keep = FALSE,
  assays.keep = "counts",
  sparsify = c("auto", "all", "none"),
  verbose = TRUE
)
```

## Arguments

- sce:

  SingleCellExperiment object

- rowData.keep:

  Either TRUE (keep all rowData columns), FALSE (remove all rowData
  columns) or a character vector with the names of the rowData columns
  to keep

- colData.keep:

  Either TRUE (keep all colData columns), FALSE (remove all colData
  columns) or a character vector with the names of the colData columns
  to keep

- metadata.keep:

  Either TRUE (keep all metadata), FALSE (remove all metadata) or a
  character vector with the names of the metadata items to keep

- assays.keep:

  Either TRUE (keep all assays), FALSE (remove all assays) or a
  character vector with the names of the assays to keep

- sparsify:

  Whether to convert assay matrices to sparse format. Either "all",
  "none" or "auto" (default) to only convert those matrices that will
  result in a size reduction

- verbose:

  Whether to print status messages

## Value

SingleCellExperiment object

## Examples

``` r

sce <- splatSimulate(verbose = FALSE)
sce.min <- minimiseSCE(sce, verbose = FALSE)
object.size(sce)
#> 46071400 bytes
object.size(sce.min)
#> 5509392 bytes
```
