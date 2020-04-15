#' Predict Method for flac Fits
#'
#' Obtains predictions from a fitted \code{flac} object. 
#' 
#' If \code{newdata} is omitted the predictions are based on the data used for the fit. 
#' 
#' 
#' @param lfobject A fitted object of class \code{flac}.
#'
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. 
#' If omitted, the fitted linear predictors are used.  
#' @param type The type of prediction required. The default is on the scale of the linear predictors. 
#' The alternative \code{"response} gives the predicted probabilities. 
#'
#' @return A vector or matrix of predictions.
#'
#' @rdname predict.flac
#' @exportS3Method predict flac
predict.flac <- function (flicobject, newdata, type = c("link", "response")) {
  type <- match.arg(type)
  names(flicobject)[3] <- "linear.predictors"
  names(flicobject)[2] <- "predict"
  predict.logistf(flicobject, newdata, type,flic=FALSE)
} 