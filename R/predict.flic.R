#' @exportS3Method predict flic
predict.flic <- function (flicobject, newdata, type = c("link", "response")) {
  type <- match.arg(type)
  names(flicobject)[3] <- "flic.linear.predictors"
  names(flicobject)[2] <- "flic.predict"
  flicobject$flic.coefficients <- flicobject$coefficients
  predict.logistf(flicobject, newdata, type,flic=TRUE)
} 