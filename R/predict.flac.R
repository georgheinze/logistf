#' @exportS3Method predict flac
predict.flac <- function (flicobject, newdata, type = c("link", "response")) {
  type <- match.arg(type)
  names(flicobject)[3] <- "linear.predictors"
  names(flicobject)[2] <- "predict"
  predict.logistf(flicobject, newdata, type,flic=FALSE)
} 