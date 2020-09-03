#' Predict Method for flic Fits
#'
#' Obtains predictions from a fitted \code{flic} object. 
#' 
#' If \code{newdata} is omitted the predictions are based on the data used for the fit. 
#' 
#' @param object A fitted object of class \code{flic}.
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.  
#' @param type The type of prediction required. The default is on the scale of the linear predictors. The alternative \code{response} gives the predicted probabilities. 
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A vector or matrix of predictions
#' @rdname predict.flic
#' @exportS3Method predict flic
predict.flic <- function(object, newdata, type = c("link", "response"), ...){
  type <- match.arg(type)
  names(object)[3] <- "flic.linear.predictors"
  names(object)[2] <- "flic.predict"
  object$flic.coefficients <- object$coefficients
  predict.logistf(object, newdata, type,flic=TRUE)
} 