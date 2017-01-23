# 0.99.8

* Fix bug that meant non-linear path factors weren't stored in output

# 0.99.7

* Small changes to avoid NOTEs and reduce check time

# 0.99.6

* Add installation to vignette
* Add detail about counts matrix to vignette
* Replace 1:x with seq_len/seq_along

# 0.99.5

* Set R_TESTS environment

# 0.99.4

* Version bump to start build

# splatter 0.99.3

* Fix to match Bioconductor version of scDD
* Add logo to repository

# splatter 0.99.2

* Add rownames, colnames to matrices in splatSimulate, lunSimulate
* Bump R version to 3.4

# splatter 0.99.1

* Address Biocondutor build warnings, notes

# splatter 0.99.0

* Submit to Bioconductor

# splatter 0.12.1

* Fix bug in compareSCESets
* Dataset order is now kept in plots

# splatter 0.12.0

* Add compareSCESets function
* Update vignette
* Fix LunParams validity bug
* Add logo

# splatter 0.11.1

* Fix bug in splatSimulatePaths that caused NAs

# splatter 0.11.0

* Make compatible with Bioconductor 3.4

# splatter 0.10.1

* Fix error for SCESets in lun2Estimate, scDDEstimate

# splatter 0.10.0

* Add listSims function
* Add vignette

# splatter 0.9.0

* Add scDD simulation

# splatter 0.8.0

* Add Lun2 simulation

# splatter 0.7.0

* Redesign how parameters are stored
* Each simulation now has it's own S4 Params class
* Modify exisiting simulations to use new parameter objects

# splatter 0.6.0

* Add Lun simulation
* Modify splatParams to take Lun parameters

# splatter 0.5.0

* Add simple simulation

# splatter 0.4.0

* Add splatter simulations
* Modify some parts of splatParams and fix bugs

# splatter 0.3.0

* Added parameter estimation functions

# splatter 0.2.0

* Added splatParams object
* Added functions for interacting with splatParams
