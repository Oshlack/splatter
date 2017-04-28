# Splatter

[![Travis-CI Build Status](https://travis-ci.org/Oshlack/splatter.svg?branch=master)](https://travis-ci.org/Oshlack/splatter)
[![Coverage Status](https://img.shields.io/codecov/c/github/Oshlack/splatter/master.svg)](https://codecov.io/github/Oshlack/splatter?branch=master)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Oshlack/splatter?branch=master&svg=true)](https://ci.appveyor.com/project/Oshlack/splatter)

![Splatter logo](https://s16.postimg.org/xc6u52b0l/splatter_logo_small.png)

Splatter is an R package for the simple simulation of single-cell RNA sequencing
data. Splatter provides a common interface to multiple simulations that have:

* Functions for estimating simulation parameters
* Objects for storing those parameters
* Functions for simulating counts using those parameters

Splatter is built on top of [`scater`][scater] and stores simulations in
`SCESet` objects. Splatter also has functions for comparing simulations and real
datasets.

## Installation.

### Development version

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

## Getting started

Once installed the best place to get started is the vignette. For most users
the most convient way to access this is online [here][vignette].

Alternatively, if you chose to build the vignette, you can load Splatter, then
browse the vignettes:

```{r}
library(splatter)
browseVignettes("splatter")
```

This is a detailed document that introduces the main features of Splatter.

[scater]: https://github.com/davismcc/scater
[contrib]: https://github.com/Bioconductor/Contributions/issues/209
[bioc]: https://bioconductor.org/packages/devel/bioc/html/splatter.html
[vignette]: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html
