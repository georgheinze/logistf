#' @exportS3Method summary flac
summary.flac <- function(object, ...)
{
  cat("Firth's logistic regression with added covariate\n\n")
  summary.logistf(object)
}