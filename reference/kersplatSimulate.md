# Kersplat simulation

Simulate scRNA-seq count data using the Kersplat model

## Usage

``` r
kersplatSimulate(
  params = newKersplatParams(),
  sparsify = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- params:

  KersplatParams object containing simulation parameters.

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment containing simulated counts and intermediate values

## Details

This functions is for simulating data in a single step. It consists of a
call to
[`kersplatSetup`](http://oshlacklab.com/splatter/reference/kersplatSetup.md)
followed by a call to
[`kersplatSample`](http://oshlacklab.com/splatter/reference/kersplatSample.md).
Please see the documentation for those functions for more details of the
individual steps.

## See also

[`kersplatSetup`](http://oshlacklab.com/splatter/reference/kersplatSetup.md),
[`kersplatSample`](http://oshlacklab.com/splatter/reference/kersplatSample.md)

## Examples

``` r

if (requireNamespace("igraph", quietly = TRUE)) {
    sim <- kersplatSimulate
}
```
