## Version 1.7.1 (2018-11-17)

* Fix bugs in vignette

## Version 1.7.0 (2018-11-01)

Bioconductor 3.9 devel

# Version 1.6.0 (2018-11-01)

Bioconductor 3.8 release

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
* Change varible name in vignette for compatibility with scater
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

* Modify how Lun2Params stores gene paramters to use data.frame
* Move sampling of genes/cells to lun2Simulate
* Return to old Lun2 nGenes estimate

## Version 0.99.12 (2017-03-22)

* Add diffSCESets function
* Update compareSCESets plots
* Modify Lun2 nGenes estimate
* Modify how addFeatureStats names columns
* Add infinte bcv.df warning to splatSimulate

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
* Modify exisiting simulations to use new parameter objects

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
