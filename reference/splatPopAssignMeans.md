# Sample expression mean and variance for each gene

A mean and coefficient of variation is assigned to each gene by sampling
from gamma distributions parameterized from real data in
`splatPopEstimate`. The cv gamma distributions are binned by gene mean
because the distribution of variance in real data is not independent
from the mean. The degree of similarity between individuals can be
further tuned using the similarity.scale parameter in `SplatPopParams`.

## Usage

``` r
splatPopAssignMeans(params, key)
```

## Arguments

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- key:

  Partial splatPop key data.frame.

## Value

The key updated with assigned means and variances.
