# logistf 1.26.1

* LINPACK routines replaced by LAPACK routines


# logistf 1.26.0

* `forward()` and `backward()` now require the dataset as an argument.
* Using `forward()` and `backward()` to select among factor variables is no longer supported.

# logistf 1.25.0

* Fixed memory bug with large data sets in `flic()` (#45).
* Added support for `emmeans`: `emmeans::emmeans()` now works with `logistf` objects (#51).


# logistf 1.24.1

* Added a `NEWS.md` file to track changes to the package.
* Added `terms()` to extract terms objects from logistf, flic and flac objects (#32).
* `logistf()` now issues warnings in case of non-convergence. 
* `logistf()` more carefully checks its response argument to ensure that they're specified correctly.
* `logistf()` has gained one new control parameter, `modcontrol`, that controls additional parameters for fitting, e.g. degree of penalization in the fitting process. Default is `logistf.mod.control()`.
* Fixed a bug in `backward()` when `data` is called data (#34) and when function is called with factors as covariates (#34).
* Updated `predict()` function for easier handling.
* `logistf()`, `flic()` and `flac()` are now compatible with `termplot()` because of the added prediction type `terms`
in `predict.logistf()`, `predict.flic()` and `predict.flac()`.
* Added a new fitting method: Iteratively reweighted least squares (IRLS). 
* Added `na.action` option to `logistf()` and `predict.logistf()`.
* Updated `flic()` and `flac()` to be called with interactions in formula.
