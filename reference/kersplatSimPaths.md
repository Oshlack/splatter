# Simulate Kersplat paths

Simulate gene means for each step along each path of a Kersplat
simulation

## Usage

``` r
kersplatSimPaths(params, verbose)
```

## Arguments

- params:

  KersplatParams object containing simulation parameters.

- verbose:

  logical. Whether to print progress messages

## Value

KersplatParams object with path means

## Details

The method of simulating paths is inspired by the method used in the
PROSSTT simulation. Changes in expression are controlled by
`paths.nPrograms` regulatory programs. Each of the regulatory genes in
the gene network has some association with each program. This is
analogous to there being changes in the environment (the programs) which
are sensed by receptors (regulatory genes) and cause changes in
expression downstream. For each path a random walk is generated for each
program and the changes passed on to the regulatory genes. At each step
the changes propagate through the network according to the weights on
edges between genes. This algorithm is fairly simple but should result
in correlation relationships between genes. However it is likely to be
improved and adjusted in the future.

The path structure itself is specified by the `paths.design` parameter.
This is a `data.frame` with three columns: "Path", "From", and "Steps".
The Path field is an ID for each path while the Steps field controls the
length of each path. Increasing the number of steps will increase the
difference in expression between the ends of the paths. The From field
sets the originating point of each path. For example a From of `0, 0, 0`
would indicate three paths from the origin while a From of `0, 1, 1`
would give a branching structure with Path 1 beginning at the origin and
Path 2 and Path 3 beginning at the end of Path 1.

## References

Papadopoulos N, Parra RG, Söding J. PROSSTT: probabilistic simulation of
single-cell RNA-seq data for complex differentiation processes.
Bioinformatics (2019). <https://doi.org/10.1093/bioinformatics/btz078>.
