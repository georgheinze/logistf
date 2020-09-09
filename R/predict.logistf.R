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
#' The alternative \code{response} gives the predicted probabilities. 
#' @param flic If \code{TRUE}(default = \code{FALSE}), predictions are computed with intercept correction.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A vector or matrix of predictions.
#'
#' @rdname predict.logistf
#' @exportS3Method predict logistf
predict.logistf <- function (object, newdata, type = c("link", "response"), flic=FALSE, ...) 
{
  type <- match.arg(type)
  if (missing(newdata)) {#no data - return linear.predictors or response according to type
    if (flic) {
      #check if flic=TRUE was set in object
      if(object$flic) {
        pred <- switch(type, link = object$linear.predictors, response = object$predict)
      }
      #if intercept is not already altered refit the model:
      else {
        message("predict called with flic=TRUE but logistf-object was called with flic=FALSE: refitting model for predictions")
        object.flic <- update(object, flic=TRUE)
        pred <- switch(type, link = object.flic$linear.predictors, response = object.flic$predict)
      }
    }
    else {
      pred <- switch(type, link = object$linear.predictors, response = object$predict)
    }
  }
  else {
    newlin <- c(object$coefficients[-1] %*% t(newdata))
    if (flic) {
      if(object$flic) {
        pred <- switch(type, link = object$coefficients[1]+newlin, 
                   response = 1/(1+exp(-(object$coefficients[1]+newlin))) )
      }
      else{
        message("predict called with flic=TRUE but logistf-object was called with flic=FALSE: refitting model for predictions")
        object.flic <- update(object, flic=TRUE)
        pred <- switch(type, link = object.flic$coefficients[1]+newlin, 
                       response = 1/(1+exp(-(object.flic$coefficients[1]+newlin))))
      }
    
    }
    else {
      pred <- switch(type, link = object$coefficients[1]+newlin , 
                     response = 1/(1+exp(-(object$coefficients[1]+newlin))))
    }
  }
  pred
}