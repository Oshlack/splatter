---
title: "Splatter"
---

<!-- badges: start -->
[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Codecov test coverage](https://codecov.io/gh/Oshlack/splatter/graph/badge.svg)](https://app.codecov.io/gh/Oshlack/splatter)
[![R-CMD-check-bioc](https://github.com/Oshlack/splatter/actions/workflows/check.yml/badge.svg)](https://github.com/Oshlack/splatter/actions/workflows/check.yml)
[![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/splatter.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/splatter)
[![Bioc devel status](http://www.bioconductor.org/shields/build/devel/bioc/splatter.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/splatter)
[![Bioc downloads rank](https://bioconductor.org/shields/downloads/release/splatter.svg)](http://bioconductor.org/packages/stats/bioc/splatter/)
[![Bioc support](https://bioconductor.org/shields/posts/splatter.svg)](https://support.bioconductor.org/tag/splatter)
[![Bioc history](https://bioconductor.org/shields/years-in-bioc/splatter.svg)](https://bioconductor.org/packages/release/bioc/html/splatter.html#since)
[![Bioc last commit](https://bioconductor.org/shields/lastcommit/devel/bioc/splatter.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/splatter/)
[![Bioc dependencies](https://bioconductor.org/shields/dependencies/release/splatter.svg)](https://bioconductor.org/packages/release/bioc/html/splatter.html#since)
<!-- badges: end -->

![Splatter logo](vignettes/splatter-logo-small.png)

Splatter is an R package for the simple simulation of single-cell RNA sequencing data.
Splatter provides a common interface to multiple simulations that have:

* Functions for estimating simulation parameters
* Objects for storing those parameters
* Functions for simulating counts using those parameters

Splatter is built on top of several [Bioconductor](bioc-home) packages and stores simulations in [`SingleCellExperiment`][SCE] objects.
Splatter also has functions for comparing simulations and real datasets.

## Installation.

Splatter is available from [Bioconductor][bioc] for R >=3.4.

It can be installed from Bioconductor with:

```r
if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("splatter")
```

If you wish to build a local version of the vignette use:

```r
BiocManager::install("splatter", build_vignettes=TRUE)
```

This will also build the vignette and install all suggested dependencies (which aren't required for core functionality).

## Getting started

Once installed the best place to get started is the vignette.
For most users the most convenient way to access this is online [here][vignette].
To get  started with population scale simulations, see the splatPop vignette [here][splatpopvignette].

Alternatively, if you chose to build the vignette, you can load Splatter, then browse the vignettes:

```r
library(splatter)
browseVignettes("splatter")
```

This is a detailed document that introduces the main features of Splatter.

## Citing Splatter

If you use Splatter please cite our paper ["Zappia L, Phipson B, Oshlack A. Splatter: Simulation Of Single-Cell RNA Sequencing Data. Genome Biology. 2017; doi:10.1186/s13059-017-1305-0"][paper].

```bibtex
@Article{,
  author = {Luke Zappia and Belinda Phipson and Alicia Oshlack},
  title = {Splatter: simulation of single-cell RNA sequencing data},
  journal = {Genome Biology},
  year = {2017},
  url = {http://dx.doi.org/10.1186/s13059-017-1305-0},
  doi = {10.1186/s13059-017-1305-0},
}
```

If you use the splatPop functions, please also cite ["Azodi CB, Zappia L, Oshlack  A, McCarthy DJ. splatPop: simulating population scale single-cell RNA sequencing data. Genome Biology. 2021; doi:10.1186/s13059-021-02546-1"][splatpoppaper].

```bibtex
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
[bioc-home]: https://www.bioconductor.org/
[vignette]: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html
[splatpopvignette]: http://www.bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatPop.html
[paper]: http://dx.doi.org/10.1186/s13059-017-1305-0
[splatpoppaper]: http://dx.doi.org/10.1186/s13059-021-02546-1
