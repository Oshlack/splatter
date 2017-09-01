[![Travis-CI Build Status](https://travis-ci.org/Oshlack/splatter.svg?branch=master)](https://travis-ci.org/Oshlack/splatter)
[![Coverage Status](https://img.shields.io/codecov/c/github/Oshlack/splatter/master.svg)](https://codecov.io/github/Oshlack/splatter?branch=master)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Oshlack/splatter?branch=master&svg=true)](https://ci.appveyor.com/project/Oshlack/splatter)
[![Bioconductor History](https://bioconductor.org/shields/years-in-bioc/splatter.svg)](https://bioconductor.org/packages/devel/bioc/html/splatter.html)
[![Bioconductor Status](https://bioconductor.org/shields/build/devel/bioc/splatter.svg)](https://bioconductor.org/packages/devel/bioc/html/splatter.html)
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/splatter.svg)](https://bioconductor.org/packages/devel/bioc/html/splatter.html)

![Splatter logo](https://s16.postimg.org/xc6u52b0l/splatter_logo_small.png)

Splatter is an R package for the simple simulation of single-cell RNA sequencing
data. Splatter provides a common interface to multiple simulations that have:

* Functions for estimating simulation parameters
* Objects for storing those parameters
* Functions for simulating counts using those parameters
* Functions for comparing simulations and real datasets

Splatter is built on top of [`scater`][scater] and stores simulations in
`SCESet` objects. Splatter also has functions for comparing simulations and real
datasets.

**!Please note!** This site provides documentation for the development version
of Splatter. For details on the current release please refer to
https://bioconductor.org/packages/splatter.

## Installation.

### Release version

Splatter has been accepted into the latest version of [Bioconductor][bioc]
and hence requires the latest version of R (>=3.4).

If you have these installed Splatter can be installed from Bioconductor using
`biocLite`:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("splatter")
```

If you wish to build a local version of the vignette use:

```{r}
biocLite("splatter", build_vignettes=TRUE)
```

This will also build the vignette and install all suggested dependencies (which
aren't required for core functionality). Building the vignette may sometimes 
fail when run from the command line, if this happens try running the install
command in RStudio.

### Development version

If you want to try the [development version][deve] this can also be installed
from Bioconductor:

```{r}
library(BiocInstaller)
useDevel()
biocValid()              # checks for out of date packages
biocLite()               # (optional) updates out of date packages
biocLite("splatter")
```

Depending on the current release cycle you may also need to install the
development version of R. See [here][bioc-install] for more details.

Alternatively the development version can be installed directly from
[Github][github]:

```{r}
devtools::install("Oshlack/splatter")
```

[scater]: https://github.com/davismcc/scater
[bioc]: https://bioconductor.org/packages/splatter
[devel]: https://bioconductor.org/packages/devel/bioc/html/splatter.html
[bioc-install]: https://www.bioconductor.org/developers/how-to/useDevel/
[github]: https://github.com/Oshlack/splatter
