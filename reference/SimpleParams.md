# The SimpleParams class

S4 class that holds parameters for the simple simulation.

## Parameters

The simple simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate.

- `nCells`:

  The number of cells to simulate.

- `[seed]`:

  Seed to use for generating random numbers.

- `mean.shape`:

  The shape parameter for the mean gamma distribution.

- `mean.rate`:

  The rate parameter for the mean gamma distribution.

- `[count.disp]`:

  The dispersion parameter for the counts negative binomial
  distribution.

The parameters not shown in brackets can be estimated from real data
using
[`simpleEstimate`](http://oshlacklab.com/splatter/reference/simpleEstimate.md).
For details of the simple simulation see
[`simpleSimulate`](http://oshlacklab.com/splatter/reference/simpleSimulate.md).
