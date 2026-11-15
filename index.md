![Splatter logo](articles/splatter-logo-small.png)

Splatter logo

Splatter is an R package for the simple simulation of single-cell RNA
sequencing data. Splatter provides a common interface to multiple
simulations that have:

- Functions for estimating simulation parameters
- Objects for storing those parameters
- Functions for simulating counts using those parameters

Splatter is built on top of several
[Bioconductor](https://www.bioconductor.org/) packages and stores
simulations in
[`SingleCellExperiment`](https://github.com/drisso/SingleCellExperiment)
objects. Splatter also has functions for comparing simulations and real
datasets.

**!Please note!** This site provides documentation for the development
version of Splatter. For details on the current release please refer to
the [Bioconductor package page](http://oshlacklab.com/splatter/bioc).

## Installation.

Splatter is available from
[Bioconductor](https://bioconductor.org/packages/devel/bioc/html/splatter.html)
for R \>=3.4.

It can be installed from Bioconductor with:

\
`if`` ``(``!`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"BiocManager"``, quietly``=``TRUE``)``)`` ``{`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"BiocManager"``)`\
`}`\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"splatter"``)`

If you wish to build a local version of the vignette use:

\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"splatter"``, build_vignettes``=``TRUE``)`

This will also build the vignette and install all suggested dependencies
(which aren’t required for core functionality).

## Getting started

Once installed the best place to get started is the vignette. For most
users the most convenient way to access this is online
[here](https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html).
To get started with population scale simulations, see the splatPop
vignette
[here](http://www.bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatPop.md).

Alternatively, if you chose to build the vignette, you can load
Splatter, then browse the vignettes:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`splatter`](https://bioconductor.org/packages/splatter/)`)`\
[`browseVignettes`](https://rdrr.io/r/utils/browseVignettes.html)`(``"splatter"``)`

This is a detailed document that introduces the main features of
Splatter.

## Citing Splatter

If you use Splatter please cite our paper [“Zappia L, Phipson B, Oshlack
A. Splatter: Simulation Of Single-Cell RNA Sequencing Data. Genome
Biology. 2017;
doi:10.1186/s13059-017-1305-0”](http://dx.doi.org/10.1186/s13059-017-1305-0).

``` bibtex
@Article{,
  author = {Luke Zappia and Belinda Phipson and Alicia Oshlack},
  title = {Splatter: simulation of single-cell RNA sequencing data},
  journal = {Genome Biology},
  year = {2017},
  url = {http://dx.doi.org/10.1186/s13059-017-1305-0},
  doi = {10.1186/s13059-017-1305-0},
}
```

If you use the splatPop functions, please also cite [“Azodi CB, Zappia
L, Oshlack A, McCarthy DJ. splatPop: simulating population scale
single-cell RNA sequencing data. Genome Biology. 2021;
doi:10.1186/s13059-021-02546-1”](http://dx.doi.org/10.1186/s13059-021-02546-1).

``` bibtex
@Article{,
  author = {Christina B Azodi and Luke Zappia and Alicia Oshlack and Davis J McCarthy},
  title = {splatPop: simulating population scale single-cell RNA sequencing data},
  journal = {Genome Biology},
  year = {2021},
  url = {http://dx.doi.org/10.1186/s13059-021-02546-1},
  doi = {10.1186/s13059-021-02546-1},
}
```
