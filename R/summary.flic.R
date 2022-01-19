 #' @exportS3Method summary flic
summary.flic <- function(object, ...)
{ 
  cat("Firth's logistic regression with intercept correction\n\n")
  summary.logistf(object)
}