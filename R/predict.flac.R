#' Predict Method for flac Fits
#'
#' Obtains predictions from a fitted \code{flac} object. 
#' 
#' If \code{newdata} is omitted the predictions are based on the data used for the fit. 
#' 
#' @param object A fitted object of class \code{flac}.
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.  
#' @param type The type of prediction required. The default is on the scale of the linear predictors. 
#' The alternative \code{response} gives the predicted probabilities. Type \code{terms} returns a matrix with the fitted
#' values of each term in the formula on the linear predictor scale.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A vector or matrix of predictions.
#' @rdname predict.flac
#' @exportS3Method predict flac
predict.flac <- function(object, newdata, type = c("link", "response", "terms"), ...){
  type <- match.arg(type)
  object$flic <- FALSE
  predict.logistf(object, newdata, type,flic=FALSE)
} 