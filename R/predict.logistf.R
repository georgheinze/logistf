#' Predict Method for logistf Fits
#'
#' Obtains predictions from a fitted \code{logistf} object. 
#' 
#' If \code{newdata} is omitted the predictions are based on the data used for the fit. 
#' 
#' 
#' @param object A fitted object of class \code{logistf}.
#'
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. 
#' If omitted, the fitted linear predictors are used.  
#' @param type The type of prediction required. The default is on the scale of the linear predictors. 
#' The alternative \code{response} gives the predicted probabilities. Type \code{terms} returns a matrix with the fitted
#' values of each term in the formula on the linear predictor scale.
#' @param flic If \code{TRUE}(default = \code{FALSE}), predictions are computed with intercept correction.
#' @param se.fit  If \code{TRUE}(default = \code{FALSE}) standard errors are computed.
#' @param reference  A named vector of reference values for each variable for \code{type="terms"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A vector or matrix of predictions.
#'
#' @rdname predict.logistf
#' @exportS3Method predict logistf

predict.logistf <- function (object, newdata, type = c("link", "response", "terms"), flic=FALSE, se.fit = FALSE, reference,...) 
{
  predict_terms <- function(object, model_matrix = NULL){
    if(is.null(model_matrix)){
      mm <- model.matrix(object$formula, object$model)
    } else {
      mm <- model_matrix
    }
      aa <- attr(mm, "assign")
      ll <- attr(terms(object), "term.labels")
      ll <- c("(Intercept)", ll)
      aaa <- factor(aa, labels = ll)
      asgn <- split(order(aa), aaa)
      asgn$"(Intercept)" <- NULL
      avx <- colMeans(mm)
      beta <- object$coefficients
      termsconst <- sum(avx * beta)
      nterms <- length(asgn)
      predictor <- matrix(ncol = nterms, nrow = NROW(mm))
      dimnames(predictor) <- list(rownames(mm), names(asgn))
      X <- sweep(mm, 2L, avx, check.margin = FALSE)
      for (i in seq.int(1L, nterms, length.out = nterms)) {
          iipiv <- asgn[[i]]
          predictor[, i] <- X[, iipiv, drop = FALSE] %*% beta[iipiv]
      }
      if(se.fit){
        se <- matrix(ncol = nterms, nrow = NROW(mm))
        dimnames(se) <- list(rownames(mm), names(asgn))
        Terms <- delete.response(terms(object))
        m <- model.frame(Terms, data.frame(t(reference)))
        reference <- model.matrix(Terms, m)
        for(t in attr(terms(object), "term.labels")){
          ind <- asgn[[t]]
          diffs <- mm[,ind,drop=FALSE]-reference[,ind]
          se_t <- apply(diffs, 1, function(x){
            t(x) %*% object$var[ind,ind] %*% x
          })
          se_t <- sqrt(se_t)
          se[, t] <- se_t
        }
      }
      attr(predictor, "constant") <- termsconst
      if (se.fit) {
        return(list(predictor = predictor, se.fit = se))
      } else {
        return(predictor)
      }
  }
  type <- match.arg(type)
  X <- model.matrix(object$formula, object$model)
  if(type == "terms" && missing(reference)){
    stop("Please provide a named vector of reference values for each variable for type=terms.")
  }
  if (missing(newdata)) {#no data - return linear.predictors or response according to type
    if (flic) {
      #check if flic=TRUE was set in object
      if(!object$flic) {
        message("predict called with flic=TRUE but logistf-object was called with flic=FALSE: refitting model for predictions")
        object <- update(object, flic=TRUE)
      }
    }
    pred <- switch(type, link = object$linear.predictors, response = object$predict, terms = predict_terms(object))
    if(se.fit && type!="terms"){
      se <- apply(X, 1, function(x){
              t(x) %*% object$var %*% x
      })
      se <- sqrt(se)
      se <- unname(se)
      if(type == "response"){
        ci_lower <- pred - 1.96*se
        ci_upper <- pred + 1.96*se
      }
    }
  }
  else {
    Terms <- delete.response(terms(object))
    m <- model.frame(Terms, newdata)
    if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl,m)
    X <-  model.matrix(Terms, m)
    if (flic) {
      if(!object$flic) {
        message("predict called with flic=TRUE but logistf-object was called with flic=FALSE: refitting model for predictions")
        object <- update(object, flic=TRUE, pl=FALSE)
      }
    }
    pred <- switch(type, link = (object$coefficients %*% t(X)) , 
                     response = (1/(1+exp(-(object$coefficients %*% t(X))))), 
                     terms = predict_terms(object))
    if(se.fit && type!="terms"){
      se <- apply(X, 1, function(x){
              t(x) %*% object$var %*% x
      })
      se <- sqrt(se)
      se <- unname(se)
      if(type == "response"){
        ci_lower <- pred - 1.96*se
        ci_upper <- pred + 1.96*se
      }
    }
  }
  if (se.fit) {
    if (type=="terms"){
      names(pred$predictor) <- NULL
      names(pred$se.fit) <- NULL
      return(list(fit = pred$predictor, se.fit = pred$se.fit))
    } else if(type=="response"){
      return(list(fit = (pred), lower = (ci_lower), upper=(ci_upper)))
    } else {
      return(list(fit = (pred), se.fit = (se)))
    }
  } else {
    return((pred))
  }
}