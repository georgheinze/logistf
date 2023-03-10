
#' Recover data method for logistf objects
#' 
#' @description This function provides a \code{recover_data} method for \code{logistf} objects. This is required for \code{emmeans} support.
#' 
#' @param object a \code{logistf} object.
#' @param frame the model frame.
#' @param ... additional arguments.
#'
#' @export recover_data.logistf

recover_data.logistf <- function(object, frame = object$model, ...){
    fcall = object$call
    emmeans::recover_data(fcall, delete.response(terms(object)),
                          object$na.action, frame = frame, ...)
}


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
    
    misc = emmeans::.std.link.labels(list(link = "logit", family = "binomial"), list())
    
    list(X = X, bhat = bhat, nbasis = nbasis, V = V,
         dffun = dffun, dfargs = dfargs, misc = misc)
}


