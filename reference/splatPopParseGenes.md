# Generate population key matrix from random or gff provided gene information

Generate population key matrix from random or gff provided gene
information

## Usage

``` r
splatPopParseGenes(params, gff)
```

## Arguments

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- gff:

  Either NULL or a data.frame object containing a GFF/GTF file.

## Value

The Partial splatPop key data.frame.
