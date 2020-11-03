# logistf

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
needed to arrive at the maximum and much more. Furthermore, specific methods for the resulting object are supplied. 
Additionally, a function to plot profiles of the penalized likelihood function and a function to perform penalized 
likelihood ratio tests have been included.

```r
data(sex2)
lf <- logistf(formula = case ~ age + oc + vic + vicl + vis + dia, data = sex2)
summary(lf)
```

<!-- badges: start -->
[![R build status](https://github.com/georgheinze/logistf/workflows/R-CMD-check/badge.svg)](https://github.com/georgheinze/logistf/actions)
<!-- badges: end -->
