# logistf 1.24.1

* Added a `NEWS.md` file to track changes to the package.
* Added `terms()` to extract terms objects from logistf, flic and flac objects (#32).
* `logistf()` now issues warnings in case of non-convergence. 
* `logistf()` more carefully checks its response argument to ensure that they're specified correctly.
* `logistf()` has gained one new parameter, `tau`, that controls the degree of penalisation in the fitting process.
* Fixed a bug in `backward()` when `data` is called data (#34).
* Updated `predict()` function for easier handling.
* Fixed a bug in `backward()` when function is called with factors as covariates (#34).
* `logistf()`, `flic()` and `flac()` are now compatible with `termplot()` because of the added prediction type `terms`
in `predict.logistf()`, `predict.flic()` and `predict.flac()`.
* Added a new fitting method: Iteratively reweighted least squares (IRLS). 
* Added `na.action` option to `logistf()` and `predict.logistf()`.
* Moved `terms.fit` and  `tau` parameters of `logistf()` to `logistf.fit.control()`.
* Added simple tests for `logistf`, `flic` and `flac`.
