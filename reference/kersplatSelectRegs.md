# Select Kersplat regulators

Select regulator genes in the gene network for a Kersplat simulation

## Usage

``` r
kersplatSelectRegs(params, verbose)
```

## Arguments

- params:

  KersplatParams object containing simulation parameters.

- verbose:

  logical. Whether to print progress messages

## Value

KersplatParams object with gene regulators

## Details

Regulators are randomly selected, weighted according to the difference
between their out degree and in degree. This is an arbitrary weighting
and may be improved or replace in the future.
