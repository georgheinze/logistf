# logistf

<!-- badges: start -->
[![R build status](https://github.com/georgheinze/logistf/workflows/RCMDcheck/badge.svg)](https://github.com/georgheinze/logistf/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/logistf)](https://cran.r-project.org/package=logistf)
[![R-CMD-check](https://github.com/georgheinze/logistf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/georgheinze/logistf/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The package logistf provides a comprehensive tool to facilitate the application of Firth’s modified
score procedure in logistic regression analysis.

## Installation
```r
# Install logistf from CRAN
install.packages("logistf")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("georgheinze/logistf")
```

## Usage
The call of the main function of the library follows the structure of the standard functions as lm
or glm, requiring a data.frame and a formula for the model specification. The resulting object belongs to the new class logistf, 
which includes penalized maximum likelihood ('Firth-Logistic'- or 'FL'-type) logistic regression parameters, standard errors, 
confidence limits, p-values, the value of the maximized penalized log likelihood, the linear predictors, the number of iterations 
needed to arrive at the maximum and much more. Furthermore, specific methods for the resulting object are supplied. The two modifications of
FL: FLIC and FLAC have been implemented. 
A function to generate and plot profiles of the penalized likelihood function and a function to perform penalized 
likelihood ratio tests are available.

```r
data(sex2)
lf <- logistf(formula = case ~ age + oc + vic + vicl + vis + dia, data = sex2)
summary(lf)
```

## Acknowledgment

This work was supported by the Austrian Science Fund (FWF) (award I 2276).
