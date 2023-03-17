
recover_data.logistf <- function(object, frame = object$model, ...){
    fcall = object$call
    emmeans::recover_data(fcall, delete.response(terms(object)),
                          object$na.action, frame = frame, ...)
}



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


