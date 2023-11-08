# Version 1.26.0 (2023-10-25)

Bioconductor 3.18 release

## Version 1.25.1 (2023-10-11)

* Fix a bug in splatSimPathDE() where DE factors were not adjusted based on the
  path origin (path.from parameter). This affected paths where the path origin
  was not the simulation origin (path.from != 0), particularly when the path DE
  was minimal.

## Version 1.25.0 (2023-04-26)

Bioconductor 3.18 devel

# Version 1.24.0 (2023-04-26)

Bioconductor 3.17 release

## Version 1.23.4 (2023-04-12)

* Replace uses of ggplot2::aes_string() with ggplot2::aes() (Fixes #75)

## Version 1.23.3 (2023-04-10)

* Fix typo in DESCRIPTION

## Version 1.23.2 (2023-04-05)

* Replace set.seed() with withr::with_seed()
* Switch scater and scuttle dependencies. scater is now suggested which saves
  a few dependencies
* Move ggplot2 to suggested for the same reason

## Version 1.23.1 (2023-01-30)

* Fix bug in splatPopSimulate() where conditional group assignments were
  incorrect when batch effects were applied (PR #161 @azodichr, fixes #158,
  fixes #160)

## Version 1.23.0 (2022-11-02)

Bioconductor 3.17 devel

# Version 1.22.0 (2022-11-02)

Bioconductor 3.16 release

## Version 1.22.1 (2023-01-30)

* Fix bug in splatPopSimulate() where conditional group assignments were
  incorrect when batch effects were applied (PR #161 @azodichr, fixes #158,
  fixes #160)

# Version 1.22.0 (2022-11-02)

Bioconductor 3.16 release

## Version 1.21.2 (2022-10-26)

* Fix a bug with non-matching rownames in splatPopSimulate() (PR #151 @azodichr,
  fixes #149)
* Fix a bug in splatPopSimulate() when batch.size is greater than the number of
  samples (PR #151 @azodichr)
* Additional splatPopSimulate tests

## Version 1.21.1 (2022-08-11)

* Fix bug in BASiCSSimulate() when spike.means is resampled (Fixes #153)

## Version 1.21.0 (2022-04-27)

Bioconductor 3.16 devel

# Version 1.20.0 (2022-04-27)

Bioconductor 3.15 release

## Version 1.19.5 (2022-04-23)

* Add starting point for dropout fit estimation (Fixes #142)
  * Idea borrowed from the modified splat implementation in the InferCNV package

## Version 1.19.4 (2022-03-11)

* New initialize method for Params objects (Fixes #132, PR #134 @wenjie2wang)

## Version 1.19.3 (2022-01-28)

* Removed duplicated code in splatSimPathCellMeans() (Fixes #131)
* Add additional fallback method for fitting dropout in splatEstimate()
  (Fixes #133)

## Version 1.19.2 (2022-01-10)

* Update the citation for the splatPop model (Fixes #128, PR #129 @azodichr)
* Update preprint DOIs for other models in listSims()
* Replace scran::computeSumFactors() with scuttle::pooledSizeFactors() in
  lun2Estimate() (Fixes #130)

## Version 1.19.1 (2021-11-02)

* Fixed bug where conditional eQTL effects were being incorrectly allocated in
  splatPop (Fixes #125, @azodichr)

## Version 1.19.0 (2021-10-27)

Bioconductor 3.15 devel

# Version 1.18.0 (2021-10-27)

Bioconductor 3.14 release

## Version 1.18.1 (2021-11-02)

* Fixed bug where conditional eQTL effects were being incorrectly allocated in
  splatPop (Fixes #125, @azodichr)

## Version 1.18.2 (2022-01-10)

* Update the citation for the splatPop model (Fixes #128, PR #129 @azodichr)
* Update preprint DOIs for other models in listSims()
* Replace scran::computeSumFactors() with scuttle::pooledSizeFactors() in
  lun2Estimate() (Fixes #130)

## Version 1.17.3 (2021-10-07)

* Improved checks for group.prob in SplatParams (Fixes #107)
* Automatically rescale group.prob during setting if it doesn't sum to 1

## Version 1.17.2 (2021-10-06)

* Fixed duplicate cell names in splatPopSimulate (Fixes #118)
* Add functionality to simulate directly from empirical values to splatPop
* Add eqtl.coreg parameter to splatPop
* Fixed a bug where too many cells were simulated in splatPop with multiple
  batches

## Version 1.17.1 (2021-05-20)

* Fix bug in validating splatPopParams

## Version 1.17.0 (2021-05-20)

Bioconductor 3.14 devel

# Version 1.16.0 (2021-05-20)

Bioconductor 3.13 release

## Version 1.16.1 (2021-05-20)

* Fix bug in validating splatPopParams

## Version 1.15.3 (2021-05-18)

* Modify sparsifyMatrices() to handle logical matrices
* Adjust checking that vectors sum to 1 (Fixes #107)
* Update NEWS.Rd for release

## Version 1.15.2 (2021-05-11)

* Updates to splatPop including conditional effects, multiplexed batch
  structure, and nCells per individual sampling. (PR #113)

## Version 1.15.1 (2020-12-01)

* Replace akima::aspline() with stats::spline(..., method = "natural") in
  buildBridges() (Fixes #108)

## Version 1.15.0 (2020-10-28)

Bioconductor 3.13 devel

# Version 1.14.0 (2020-10-28)

Bioconductor 3.12 release

## Version 1.14.1 (2020-12-01)

* Replace akima::aspline() with stats::spline(..., method = "natural") in
  buildBridges() (Fixes #108)

## Version 1.13.2 (2020-10-25)

* Add batch.rmEffect parameter to SplatParams (PR #103)
* Add minimiseSCE() function
* Automatically sparsify assays by default in all simulations
* Add checkDependencies() function
* Replace DropletUtils::downsampleMatrix() with scuttle::downsampleMatrix()

## Version 1.13.1 (2020-10-22)

* Add the splatPop simulation (PR #106)

## Version 1.13.0 (2020-04-29)

* Bioconductor 3.12 devel

# Version 1.12.0 (2020-04-29)

* Bioconductor 3.11 release

## Version 1.11.6 (2020-04-22)

* Fix warning about typo in NEWS.Rd
* Fix note about partial argument matching

## Version 1.11.5 (2020-04-20)

* Use alternative algorithm if splatEstimate dropout fitting fails (Fixes #96
  and #31)

## Version 1.11.4 (2020-04-01)

* Adjust paths example in vignette (Fixes #90)
* Use option for Kersplat warning to avoid NOTE about ':::'

## Version 1.11.3 (2020-02-20)

* Check for cycles in path.from (Fixes #91)

## Version 1.11.2 (2020-02-14)

* Update tests to work with checkmate 2.0.0
* Add suggested checks to Kersplat tests and examples

## Version 1.10.1 (2020-02-14)

* Update tests to work with checkmate 2.0.0

## Version 1.11.1 (2010-01-30)

* Replace defunct scater function in vignettes
* Update Roxygen version

## Version 1.11.0 (2019-10-30)

* Bioconductor 3.11 devel

# Version 1.10.0 (2019-10-30)

* Bioconductor 3.10 release

## Version 1.9.11 (2019-10-20)

* Update NEWS.Rd

## Version 1.9.10 (2019-10-20)

* Update scater QC function names

## Version 1.9.9 (2019-10-16)

* Add missing documentation for sce argument in getCounts
* Update maintainer email address

## Version 1.9.8 (2019-10-11)

* Add Kersplat simulation! Still experimental but is useable.
* Check for counts assay when estimating from SingleCellExperiment objects
  (Fixes #82)
* Fix where simpleSimulate stores parameters (Fixes #72)

### Version 1.9.7.9019 (2019-10-11)

* Run checks

### Version 1.9.7.9018 (2019-10-10)

* Document Splotch estimation functions
* Add Splotch to listSims
* Rename Splotch to Kersplat

### Version 1.9.7.9017 (2019-10-09)

* Document Splotch simulation functions

### Version 1.9.7.9016 (2019-10-01)

* Change the default SplotchParams cell design to place cells at end of path
* Add one time warning for Splotch simulation
* Avoid unnecessarily resetting cells.design in SplotchParams

### Version 1.9.7.9015 (2019-09-26)

* Add doublets to splotchSimulate

### Version 1.9.7.9014 (2019-09-24)

* Add ambient expression and empty cells to splotchSimulate

### Version 1.9.7.9013 (2019-09-19)

* Merge master into splotch branch

## Version 1.9.7 (2019-09-19)

* Rescale when sampling Phi in BASiCSSimulate

## Version 1.9.6 (2019-09-17)

* Adjust BASiCSSimulate to match development version of BASiCS

### Version 1.9.5.9012 (2019-09-09)

* Merge master into splotch branch

## Version 1.9.5 (2019-09-09)

* Fix bug in NAMESPACE

### Version 1.9.4.9011 (2019-08-29)

* Add BCV adjustment to Splotch simulation
* Use new exponential correction for bcv.common

### Version 1.9.4.9010 (2019-08-22)

* Add density sampling options for means and library sizes

### Version 1.9.4.9009 (2019-08-21)

* Replace library size log-normal with density and rejection sampling
* Add violins to comparison plots

### Version 1.9.4.9008 (2019-08-20)

* Merge master into splotch branch
* Fix dataset in splotchEstimate example

## Version 1.9.4 (2019-08-20)

* Fix deprecated column name in diffSCEs QQ plots
* Fix bugs where parameters were not being passed correctly in BASiCSEstimate
  and sparseDCEstimate
* Replace the sc_example_counts dataset from scater with the mockSCE function
* Tidy and improve estimation function examples

### Version 1.9.3.9008 (2019-08-20)

* Add expression outliers to Splotch

### Version 1.9.3.9007 (2019-08-14)

* Fix bug in selectFit

### Version 1.9.3.9006 (2019-08-13)

* Add splotchEstimate function

### Version 1.9.3.9005 (2019-08-08)

* Simulate counts in splotchSimulate
* Split splotchSimulate into separate functions for each stage

### Version 1.9.3.9004 (2019-08-08)

* Merge master into splotch branch

## Version 1.9.3 (2019-08-08)

* Remove deprecated scater functions

### Version 1.9.2.9003 (2019-07-17)

* Topologically sort Splotch paths
* Add library size parameters to SplotchParams
* Simulate cell means in splotchSimulate

### Version 1.9.2.9002 (2019-07-16)

* Add paths parameters to SplotchParams

### Version 1.9.2.9001 (2019-07-16)

* Add mean parameters to SplotchParams

### Version 1.9.2.9000 (2019-07-11)

* Add SplotchParams object

## Version 1.9.2 (2019-06-13)

* Add variable gene correlation plot to compareSCEs

## Version 1.9.1 (2019-06-05)

* Refactor the summariseDiff function
* Add the KS statistic to summariseDiff
* Set gam method for fitted lines in comparison plots
* Add checks for suggested packages in examples

## Version 1.9.0 (2019-05-03)

* Bioconductor 3.10 devel

# Version 1.8.0 (2019-05-03)

* Bioconductor 3.9 release

## Version 1.7.5 (2018-04-18)

* Add Splat parameters vignette
* Fix spelling

## Version 1.7.4 (2018-04-17)

* Allow SplatParams parameters to be set in any order

## Version 1.7.3 (2018-04-15)

* Minor changes to documentation
* Protect against integer overflow in simulation functions
* Fix order of groups in splatSimulate

## Version 1.7.2 (2018-12-06)

* Add ImmunoOncology biocViews tag at request of Bioconductor team

## Version 1.7.1 (2018-11-17)

* Fix bugs in vignette

## Version 1.7.0 (2018-11-01)

* Bioconductor 3.9 devel

## Version 1.6.1 (2018-12-06)

* Add ImmunoOncology biocViews tag at request of Bioconductor team

# Version 1.6.0 (2018-11-01)

* Bioconductor 3.8 release

## Version 1.5.8 (2018-10-22)

* Add WithSpikes argument when calling BASiCS::BASiCS_MCMC()

## Version 1.5.7 (2018-10-03)

* Fix BASiCSEstimate tests

## Version 1.5.4 (2018-08-30)

* Fix installation instructions

## Version 1.5.3 (2018-08-20)

* Fix bug in BASiCSEstimate when no spike-ins provided

### Version 1.5.2.9000 (2018-08-16)

* Add missing regression argument to BASiCSEstimate
* Add new BASiCS tests

## Version 1.5.2 (2018-08-16)

* Fix bug in getLNormFactors when reversing factors less than one
* Update documentation to new Roxygen version (6.1.0)
* Change variable name in vignette for compatibility with scater
* Add suggested package checks to tests

## Version 1.5.1 (2018-06-12)

* Fix normality testing error in splatEstLib
* Correct p-value cutoff in normality test
* Sample library sizes for normality testing if > 5000 cells

## Version 1.5.0 (2018-05-02)

* Bioconductor 3.7 devel

## Version 1.4.3 (2018-08-20)

* Fix bug in BASiCSEstimate when no spike-ins provided

## Version 1.4.2 (2018-08-16)

* Add missing regression argument to BASiCSEstimate
* Add new BASiCS tests

## Version 1.4.1 (2018-06-12)

* Fix normality testing error in splatEstLib
* Correct p-value cutoff in normality test
* Sample library sizes for normality testing if > 5000 cells

# Version 1.4.0 (2018-05-02)

* Bioconductor 3.7 release

## Version 1.3.6 (2018-04-30)

* Set seed in sparseDCEstimate tests

## Version 1.3.5 (2018-04-25)

* Replace dropout.present with dropout.type in SplatParams
    * Allows users to set dropout by experiment, batch, group or cell

## Version 1.3.4 (2018-04-19)

* Add option to use a normal distribution for library sizes in Splat simulations
* Add SparseDC simulation
* Rename params in metadata slot of simulation to Params for consistency
* Add checks to SplatParams to reset path.from if nGroups changes
* Improve printing of parameters stored as data.frames
* Colourise Params print output and listSims output
* Improve warning message when fitting means in splatEstimate fails
* Simplify problematic documentation links

### Version 1.3.3.9010 (2018-04-12)

* Add option to use a normal distribution for library sizes in Splat simulations

### Version 1.3.3.9000 (2018-04-12)

* Replace dropout.present with dropout.type in SplatParams
    * Allows more control over dropout.mid and dropout.shape

## Version 1.3.3 (2018-03-27)

* Fix parameter passing bug in scDDEstimate

## Version 1.3.2 (2018-01-31)

* Fix QC names that have been changed in scater
* Move scater to Imports and add scater version
* Remove lingering references to SCESets
* Add more tests

## Version 1.3.1 (2018-01-03)

* Fix error in vignette caused by changes in scater

## Version 1.3.0 (2017-10-13)

* Bioconductor 3.7 devel

## Version 1.2.2 (2018-03-27)

* Fix parameter passing bug in scDDEstimate

## Version 1.2.1 (2017-11-23)

* Fix zinbwave installation error

# Version 1.2.0 (2017-10-30)

* Bioconductor 3.6 release

## Version 1.1.8 (2017-10-13)

* Add BASiCS simulation
* Update Splatter citation
* Update Lun2 reference

## Version 1.1.7 (2017-10-05)

* Add PhenoPath simulation
* Add ZINB-WaVE simulation
* Adjust mfaSimulate output

## Version 1.1.6 (2017-10-02)

* Update scDD version
* Add mfa simulation

## Version 1.1.5 (2017-09-13)

* Convert to SingleCellExperiment

## Version 1.1.4 (2017-08-04)

* Fix group factors bug

## Version 1.1.3 (2017-07-20)

* Add verbose option to scDDEstimate
* Change "mean-dropout" to "mean-zeros" in compareSCESets

## Version 1.1.2 (2017-07-16)

* Update summariseDiff
* Update scDDEstimate, now estimates gene types
* Fix error in lun2Estimate
* Import SummarizedExperiment to avoid warnings

## Version 1.1.1 (2017-07-07)

* Add batch effects to Splat simulation
* Make Splat group assignment probabilistic
* Update SplatParams with new parameters

## Version 1.1.0 (2017-07-07)

* Bioconductor 3.6 devel

## Version 1.0.3 (2017-05-26)

* Update citation

## Version 1.0.2 (2017-05-15)

* Fix error handling when fitting means

## Version 1.0.1 (2017-04-28)

* Fix scales in some difference plots
* Fix colours in difference plots
* Fix panel legends

# Version 1.0.0 (2017-04-28)

* Bioconductor 3.5 release

## Version 0.99.15 (2017-04-14)

* Add summariseDiff function
* Add BPPARAM argument to scDDSimulate
* Adjust default Splat DE factor parameters
* Add limits to zeros diff plots
* Remove estimation of dropout.present

## Version 0.99.14 (2017-03-28)

* Add functions for making comparison panels
* Add panel section to vignette
* Change variance plot scale (for consistency)

## Version 0.99.13 (2017-03-25)

* Modify how Lun2Params stores gene parameters to use data.frame
* Move sampling of genes/cells to lun2Simulate
* Return to old Lun2 nGenes estimate

## Version 0.99.12 (2017-03-22)

* Add diffSCESets function
* Update compareSCESets plots
* Modify Lun2 nGenes estimate
* Modify how addFeatureStats names columns
* Add infinite bcv.df warning to splatSimulate

## Version 0.99.11 (2017-03-20)

* Add parallel option to lun2Estimate
* Allow non-integer library sizes in Lun2Params
* Adjust dropout eta value
* Adjust comparison plots

## Version 0.99.10 (2017-03-07)

* Improve Splat estimation procedure
* Update default Splat parameters
* Remove out.loProb Splat parameter
* Add MeanZeros plot to compareSCESets

## Version 0.99.9 (2017-02-02)

* Add addGeneLengths function
* Update scDD reference

## Version 0.99.8 (2017-01-23)

* Fix bug that meant non-linear path factors weren't stored in output

## Version 0.99.7 (2017-01-10)

* Small changes to avoid NOTEs and reduce check time

## Version 0.99.6 (2016-12-29)

* Add installation to vignette
* Add detail about counts matrix to vignette
* Replace 1:x with seq_len/seq_along

## Version 0.99.5 (2016-12-28)

* Set R_TESTS environment

## Version 0.99.4 (2016-12-23)

* Version bump to start build

## Version 0.99.3 (2016-12-21)

* Fix to match Bioconductor version of scDD
* Add logo to repository

## Version 0.99.2 (2016-12-13)

* Add rownames, colnames to matrices in splatSimulate, lunSimulate
* Bump R version to 3.4

## Version 0.99.1 (2016-12-12)

* Address Biocondutor build warnings, notes

# Version 0.99.0 (2016-12-05)

* Submit to Bioconductor

## Version 0.12.1 (2016-11-25)

* Fix bug in compareSCESets
* Dataset order is now kept in plots

## Version 0.12.0 (2016-10-25)

* Add compareSCESets function
* Update vignette
* Fix LunParams validity bug
* Add logo

## Version 0.11.1 (2016-11-23)

* Fix bug in splatSimulatePaths that caused NAs

## Version 0.11.0 (2016-11-22)

* Make compatible with Bioconductor 3.4

## Version 0.10.1 (2016-10-17)

* Fix error for SCESets in lun2Estimate, scDDEstimate

## Version 0.10.0 (2016-10-16)

* Add listSims function
* Add vignette

## Version 0.9.0 (2016-10-15)

* Add scDD simulation

## Version 0.8.0 (2016-10-15)

* Add Lun2 simulation

## Version 0.7.0 (2016-10-14)

* Redesign how parameters are stored
* Each simulation now has it's own S4 Params class
* Modify existing simulations to use new parameter objects

## Version 0.6.0 (2016-10-12)

* Add Lun simulation
* Modify splatParams to take Lun parameters

## Version 0.5.0 (2016-10-12)

* Add simple simulation

## Version 0.4.0 (2016-10-12)

* Add splatter simulations
* Modify some parts of splatParams and fix bugs

## Version 0.3.0 (2016-10-06)

* Added parameter estimation functions

## Version 0.2.0 (2016-10-06)

* Added splatParams object
* Added functions for interacting with splatParams
