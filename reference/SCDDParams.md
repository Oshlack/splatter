# The SCDDParams class

S4 class that holds parameters for the scDD simulation.

## Parameters

The SCDD simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate (not used).

- `nCells`:

  The number of cells to simulate in each condition.

- `[seed]`:

  Seed to use for generating random numbers.

- `SCdat`:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  containing real data.

- `nDE`:

  Number of DE genes to simulate.

- `nDP`:

  Number of DP genes to simulate.

- `nDM`:

  Number of DM genes to simulate.

- `nDB`:

  Number of DB genes to simulate.

- `nEE`:

  Number of EE genes to simulate.

- `nEP`:

  Number of EP genes to simulate.

- `[sd.range]`:

  Interval for fold change standard deviations.

- `[modeFC]`:

  Values for DP, DM and DB mode fold changes.

- `[varInflation]`:

  Variance inflation factors for each condition. If all equal to 1 will
  be set to `NULL` (default).

- `[condition]`:

  String giving the column that represents biological group of interest.

The parameters not shown in brackets can be estimated from real data
using
[`scDDEstimate`](http://oshlacklab.com/splatter/reference/scDDEstimate.md).
See [`simulateSet`](https://rdrr.io/pkg/scDD/man/simulateSet.html) for
more details about the parameters. For details of the Splatter
implementation of the scDD simulation see
[`scDDSimulate`](http://oshlacklab.com/splatter/reference/scDDSimulate.md).
