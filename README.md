# Splatter

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Travis-CI Build Status](https://travis-ci.org/Oshlack/splatter.svg?branch=master)](https://travis-ci.org/Oshlack/splatter)
[![Coverage Status](https://img.shields.io/codecov/c/github/Oshlack/splatter/master.svg)](https://codecov.io/github/Oshlack/splatter?branch=master)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Oshlack/splatter?branch=master&svg=true)](https://ci.appveyor.com/project/Oshlack/splatter)
![R-CMD-check](https://github.com/Oshlack/splatter/workflows/R-CMD-check/badge.svg)
[![Bioc Years](https://bioconductor.org/shields/years-in-bioc/splatter.svg)](https://bioconductor.org/packages/devel/bioc/html/splatter.html)
[![Bioc Stats](https://bioconductor.org/shields/downloads/splatter.svg)](https://bioconductor.org/packages/devel/bioc/html/splatter.html)
[![Bioc Build](https://bioconductor.org/shields/build/devel/bioc/splatter.svg)](https://bioconductor.org/packages/devel/bioc/html/splatter.html)

![Splatter logo](vignettes/splatter-logo-small.png)

Splatter is an R package for the simple simulation of single-cell RNA sequencing
data. Splatter provides a common interface to multiple simulations that have:

* Functions for estimating simulation parameters
* Objects for storing those parameters
* Functions for simulating counts using those parameters

Splatter is built on top of [`scater`][scater] and stores simulations in
[`SingleCellExperiment`][SCE] objects. Splatter also has functions for comparing
simulations and real datasets.

## Installation.

Splatter is available from [Bioconductor][bioc] for R >=3.4.

It can be installed from Bioconductor with:

```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("splatter")
```

If you wish to build a local version of the vignette use:

```{r}
BiocManager::install("splatter", build_vignettes=TRUE)
```

This will also build the vignette and install all suggested dependencies (which
aren't required for core functionality).

## Getting started

Once installed the best place to get started is the vignette. For most users
the most convenient way to access this is online [here][vignette]. To get 
started with population scale simulations, see the splatPop vignette 
[here][splatpopvignette].

Alternatively, if you chose to build the vignette, you can load Splatter, then
browse the vignettes:

```{r}
library(splatter)
browseVignettes("splatter")
```

This is a detailed document that introduces the main features of Splatter.

## Citing Splatter

If you use Splatter please cite our paper ["Zappia L, Phipson B, Oshlack A.
Splatter: Simulation Of Single-Cell RNA Sequencing Data. Genome Biology. 2017;
doi:10.1186/s13059-017-1305-0"][paper].

```
  @Article{,
    author = {Luke Zappia and Belinda Phipson and Alicia Oshlack},
    title = {Splatter: simulation of single-cell RNA sequencing data},
    journal = {Genome Biology},
    year = {2017},
    url = {http://dx.doi.org/10.1186/s13059-017-1305-0},
    doi = {10.1186/s13059-017-1305-0},
  }
```

If you use the splatPop functions, please also cite ["Azodi CB, Zappia L, Oshlack 
A, McCarthy DJ. splatPop: simulating population scale single-cell RNA sequencing 
data. Genome Biology. 2021; doi:10.1186/s13059-021-02546-1"][splatpoppaper].

```
  @Article{,
    author = {Christina B Azodi and Luke Zappia and Alicia Oshlack and Davis J McCarthy},
    title = {splatPop: simulating population scale single-cell RNA sequencing data},
    journal = {Genome Biology},
    year = {2021},
    url = {http://dx.doi.org/10.1186/s13059-021-02546-1},
    doi = {10.1186/s13059-021-02546-1},
  }
```

[scater]: https://github.com/davismcc/scater
[SCE]: https://github.com/drisso/SingleCellExperiment
[contrib]: https://github.com/Bioconductor/Contributions/issues/209
[bioc]: https://bioconductor.org/packages/devel/bioc/html/splatter.html
[vignette]: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html
[splatpopvignette]: http://www.bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatPop.html
[paper]: http://dx.doi.org/10.1186/s13059-017-1305-0
[splatpoppaper]: http://dx.doi.org/10.1186/s13059-021-02546-1
