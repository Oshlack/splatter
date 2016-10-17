# Splatter

[![Travis-CI Build Status](https://travis-ci.org/Oshlack/splatter.svg?branch=master)](https://travis-ci.org/Oshlack/splatter)
[![Coverage Status](https://img.shields.io/codecov/c/github/Oshlack/splatter/master.svg)](https://codecov.io/github/Oshlack/splatter?branch=master)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Oshlack/splatter?branch=master&svg=true)](https://ci.appveyor.com/project/Oshlack/splatter)

Splatter is an R package for the simple simulation of single-cell RNA sequencing
data. Splatter provides a common interface that has:

* Functions for estimating simulation parameters
* Objects for storing those parameters
* Functions for simulating counts using those parameters

Splatter is built on top of [`scater`][scater] and stores simulations in
`SCESet` objects.

## Installation.

Splatter can be installed from Github using `biocLite`:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("Oshlack/splatter", dependencies=TRUE, build_vignettes=TRUE)
```

This will also build the vignette and install suggested dependencies (which
aren't required for core functionality). Building the vignette may sometimes 
fail when run from the command line, if this happens try running the install
command in RStudio. If you still have problems or want a quicker install
(without the vignette) use:

```{r}
biocLite("Oshlack/splatter")
```

## Getting started

Once installed the best place to get started is the vignette. Load Splatter,
then browse the vignettes:

```{r}
library(splatter)
browseVignettes("splatter")
```
This is a detailed document that introduces the main features of Splatter.

[scater]: https://github.com/davismcc/scater