#' Predict Method for logistf Fits
#'
#' Obtains predictions from a fitted \code{logistf} object. 
#' 
#' If \code{newdata} is omitted the predictions are based on the data used for the fit. 
#' 
#' 
#' @param lfobject A fitted object of class \code{logistf}.
#'
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. 
#' If omitted, the fitted linear predictors are used.  
#' @param type The type of prediction required. The default is on the scale of the linear predictors. 
#' The alternative \code{"response} gives the predicted probabilities. 
#' @param flic If \code{TRUE}(default = \code{FALSE}), predictions are computed with intercept correction.
#'
#' @return A vector or matrix of predictions.
#'
#' @rdname predict.logistf
#' @exportS3Method predict logistf
predict.logistf <- function (lfobject, newdata, type = c("link", "response"), flic=FALSE) 
{
  type <- match.arg(type)
  if (missing(newdata)) {#no data - return linear.predictors or response according to type
    if (flic) {
      pred <- switch(type, link = lfobject$flic.linear.predictors, response = lfobject$flic.predict)
    }
    else {
      pred <- switch(type, link = lfobject$linear.predictors, response = lfobject$predict)
    }
  }
  else {
    newlin <- c(lfobject$coefficients[-1] %*% t(newdata))
    if (flic) {
    pred <- switch(type, link = lfobject$flic.coefficients[1]+newlin, 
                   response = 1/(1+exp(-(lfobject$flic.coefficients[1]+newlin))) )
    }
    else {
      pred <- switch(type, link = lfobject$coefficients[1]+newlin , 
                     response = 1/(1+exp(-(lfobject$coefficients[1]+newlin))))
    }
  }
  pred
}