#' Predict Method for flic Fits
#'
#' Obtains predictions from a fitted \code{flic} object. 
#' 
#' If \code{newdata} is omitted the predictions are based on the data used for the fit. 
#' 
#' @param object A fitted object of class \code{flic}.
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.  
#' @param type The type of prediction required. The default is on the scale of the linear predictors. 
#' The alternative \code{response} gives the predicted probabilities. Type \code{terms} returns a matrix with the fitted
#' values of each term in the formula on the linear predictor scale.
#' @param se.fit  If \code{TRUE}(default = \code{FALSE}) standard errors are computed.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A vector or matrix of predictions
#' @rdname predict.flic
#' @exportS3Method predict flic
predict.flic <- function(object, newdata, type = c("link", "response", "terms"), se.fit = FALSE, ...){
  type <- match.arg(type)
  object$flic <- TRUE
  return(predict.logistf(object, newdata, type,flic=TRUE, se.fit = se.fit, ...))
} 