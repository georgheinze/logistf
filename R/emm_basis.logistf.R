#' emm basis method for logistf objects
#' 
#' @description This function provides a \code{emm_basis} method for \code{logistf} objects. This is required for \code{emmeans} support.
#' 
#' @param object a \code{logistf} object.
#' @param trms a model terms object.
#' @param xlev a named list of character vectors giving the full set of levels to be assumed for each factor.
#' @param grid a data.frame, list or environment (or object coercible by as.data.frame to a data.frame), containing the variables in formula. Neither a matrix nor an array will be accepted.
#' @param ... additional arguments.
#' 
#' @export emm_basis.logistf

emm_basis.logistf <- function(object, trms, xlev, grid, ...) { 
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts) 
    bhat = coef(object) 
    Xmat = model.matrix(trms, data=object$model)
    V = vcov(object)
    nbasis = matrix(NA) 
    dfargs = list(df = nrow(Xmat) - ncol(Xmat))
    dffun = function(k, dfargs) dfargs$df
    list(X = X, bhat = bhat, nbasis = nbasis, V = V,
         dffun = dffun, dfargs = dfargs)
}


