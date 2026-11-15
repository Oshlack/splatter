# Get log-normal factors

Randomly generate multiplication factors from a log-normal distribution.

## Usage

``` r
getLNormFactors(n.facs, sel.prob, neg.prob, fac.loc, fac.scale)
```

## Arguments

- n.facs:

  Number of factors to generate.

- sel.prob:

  Probability that a factor will be selected to be different from 1.

- neg.prob:

  Probability that a selected factor is less than one.

- fac.loc:

  Location parameter for the log-normal distribution.

- fac.scale:

  Scale factor for the log-normal distribution.

## Value

Vector containing generated factors.
